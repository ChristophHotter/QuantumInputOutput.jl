using StaticArrays: StaticArrays, SMatrix, SVector
using FunctionWrappers: FunctionWrapper

# NOTE: FunctionWrapper is callable but NOT <: Function.
const _Callable = Union{Function,FunctionWrapper}

# ──────────────────────────────────────────────
# Dispatch helpers: symbolic vs numeric
# ──────────────────────────────────────────────

_adj(x::SQA.QField) = adjoint(x)
_adj(f::_Callable) = t -> adjoint(f(t))
_adj(x) = adjoint(x)

_post(x::BasicSymbolic) = simplify(x)
_post(x) = x

# ──────────────────────────────────────────────
# Arithmetic helpers: static vs time-dependent
# ──────────────────────────────────────────────

_is_time_dep(x) = x isa _Callable
_to_func(x) = _is_time_dep(x) ? x : (t -> x)

_isunit(x::Number) = isone(x)
_isunit(::Any) = false

_add(f::_Callable, g::_Callable) = t -> f(t) + g(t)
_add(f::_Callable, x) = iszero(x) ? f : (t -> f(t) + x)
_add(x, f::_Callable) = iszero(x) ? f : (t -> x + f(t))
_add(x, y) = iszero(x) ? y : iszero(y) ? x : x + y

_mul(f::_Callable, g::_Callable) = t -> f(t) * g(t)
_mul(s, f::_Callable) = isone(s) ? f : (t -> s * f(t))
_mul(f::_Callable, s) = isone(s) ? f : (t -> f(t) * s)
_mul(x, y) = _isunit(x) ? y : _isunit(y) ? x : x * y

# ──────────────────────────────────────────────
# FunctionWrapper with concrete return type
# ──────────────────────────────────────────────

"""
    _detect_operator_type(L, H)

Determine the concrete return type for FunctionWrapper by inspecting
the elements of L and H. Checks FunctionWrapper type parameters first
(no evaluation needed), then static element types.
Errors if only plain closures are present (the caller should thread
the type through from input SLH objects instead).
"""
function _detect_operator_type(L, H)
    # Check FunctionWrapper elements (carry explicit type info)
    for l in L
        l isa FunctionWrapper && return _fw_return_type(typeof(l))
    end
    H isa FunctionWrapper && return _fw_return_type(typeof(H))
    # Check static (non-time-dep) elements
    for l in L
        _is_time_dep(l) || return typeof(l)
    end
    _is_time_dep(H) || return typeof(H)
    # No type information available
    error(
        "Cannot determine concrete operator type: all elements are untyped closures. " *
        "Wrap at least one with FunctionWrapper{OpType, Tuple{Float64}} or include a static element.",
    )
end

_fw_return_type(::Type{FunctionWrapper{R,A}}) where {R,A} = R

function _maybe_wrap_lindblad(L::SVector{N}, ::Type{OpType}) where {N,OpType}
    if any(_is_time_dep, L)
        fw_type = FunctionWrapper{OpType,Tuple{Float64}}
        return SVector{N,fw_type}(ntuple(i -> fw_type(_to_func(L[i])), Val(N)))
    end
    return L
end

function _maybe_wrap_hamiltonian(H, has_td::Bool, ::Type{OpType}) where {OpType}
    if has_td
        return FunctionWrapper{OpType,Tuple{Float64}}(_to_func(H))
    end
    return H
end

# ──────────────────────────────────────────────
# SLH type
# ──────────────────────────────────────────────

"""
    SLH{N, ST, LT, HT}

SLH triple with scattering matrix `S`, Lindblad vector `L`, and Hamiltonian `H`.
`S` and `L` can also be vectors of scattering matrices and Lindblad terms.

See also [`▷`](@ref), [`⊞`](@ref), [`feedback`](@ref)
"""
struct SLH{N,ST,LT,HT,L}
    scattering::SMatrix{N,N,ST,L}
    lindblad::SVector{N,LT}
    hamiltonian::HT
    function SLH{N,ST,LT,HT}(
        S::SMatrix{N,N,ST,L}, lindblad::SVector{N,LT}, H::HT
    ) where {N,ST,LT,HT,L}
        return new{N,ST,LT,HT,L}(S, lindblad, H)
    end
end

"""
    _op_type(G::SLH)

Extract the operator return type from an SLH with FunctionWrapper elements.
Returns `nothing` if no FunctionWrapper type info is available.
"""
function _op_type(::SLH{N,ST,LT,HT}) where {N,ST,LT,HT}
    LT <: FunctionWrapper && return _fw_return_type(LT)
    HT <: FunctionWrapper && return _fw_return_type(HT)
    return nothing
end

# ──────────────────────────────────────────────
# Constructors
# ──────────────────────────────────────────────

# Canonical: from SMatrix + SVector (handles FunctionWrapper wrapping)
function _build_slh(S::SMatrix{N,N}, L::SVector{N}, H) where {N}
    has_td = any(_is_time_dep, L) || _is_time_dep(H)
    if has_td
        OpType = _detect_operator_type(L, H)
        return _build_slh(S, L, H, OpType)
    end
    return SLH{N,eltype(S),eltype(L),typeof(H)}(S, L, H)
end

# With explicit operator type (skips detection — used by composition operations)
function _build_slh(S::SMatrix{N,N}, L::SVector{N}, H, ::Type{OpType}) where {N,OpType}
    has_td = any(_is_time_dep, L) || _is_time_dep(H)
    if has_td
        L_w = _maybe_wrap_lindblad(L, OpType)
        H_w = _maybe_wrap_hamiltonian(H, true, OpType)
        return SLH{N,eltype(S),eltype(L_w),typeof(H_w)}(S, L_w, H_w)
    end
    return SLH{N,eltype(S),eltype(L),typeof(H)}(S, L, H)
end

# Nothing hint falls through to detection
_build_slh(S::SMatrix{N,N}, L::SVector{N}, H, ::Nothing) where {N} = _build_slh(S, L, H)

# From AbstractMatrix + AbstractVector (includes SMatrix + SVector)
function SLH(S::AbstractMatrix, L::AbstractVector, H)
    N = length(L)
    @assert size(S, 1) == N && size(S, 2) == N
    return _build_slh(SMatrix{N,N}(S), SVector{N}(L...), H)
end

# Numeric scalar S + vector L → S * I_{NxN}
function SLH(S::Number, L::AbstractVector, H)
    N = length(L)
    S_mat = SMatrix{N,N}(S * LinearAlgebra.I)
    return _build_slh(S_mat, SVector{N}(L...), H)
end

# Scalar S + scalar L → SLH{1}
function SLH(S, L, H)
    S_mat = SMatrix{1,1}(S)
    L_vec = SVector{1}(L)
    return _build_slh(S_mat, L_vec, H)
end

# Symbolic/general scalar S + vector L → S * I
function SLH(S, L::AbstractVector, H)
    N = length(L)
    S_mat = SMatrix{N,N}([i == j ? S : 0 for i = 1:N, j = 1:N])
    return _build_slh(S_mat, SVector{N}(L...), H)
end

# ──────────────────────────────────────────────
# Accessors
# ──────────────────────────────────────────────

"""
    scattering(G::SLH)

Return the scattering matrix `S` of an SLH object.
"""
scattering(G::SLH) = G.scattering

"""
    lindblad(G::SLH)

Return the Lindblad vector `L` of an SLH object.
"""
lindblad(G::SLH) = G.lindblad

"""
    hamiltonian(G::SLH)

Return the Hamiltonian `H` of an SLH object.
"""
hamiltonian(G::SLH) = G.hamiltonian


# ──────────────────────────────────────────────
# Equality
# ──────────────────────────────────────────────

function Base.isequal(G1::SLH, G2::SLH)
    isequal(scattering(G1), scattering(G2)) &&
        isequal(lindblad(G1), lindblad(G2)) &&
        isequal(hamiltonian(G1), hamiltonian(G2))
end

# ──────────────────────────────────────────────
# Matrix-vector helpers
# ──────────────────────────────────────────────

@generated function _slh_matvec(S::SMatrix{N,N}, L::SVector{N}) where {N}
    if N == 1
        return :(SVector{1}(_mul(S[1, 1], L[1])))
    end
    exprs = []
    for i = 1:N
        first_name = Symbol("tmp_$(i)_1")
        terms = [:($first_name = _mul(S[$i, 1], L[1]))]
        acc = first_name
        for j = 2:N
            tname = Symbol("tmp_$(i)_$(j)")
            push!(terms, :($tname = _add($acc, _mul(S[$i, $j], L[$j]))))
            acc = tname
        end
        push!(exprs, Expr(:block, terms..., acc))
    end
    return :(SVector($(exprs...)))
end

@generated function _slh_dot(L1::SVector{N}, L2::SVector{N}) where {N}
    if N == 1
        return :(_mul(L1[1], L2[1]))
    end
    expr = :(_mul(L1[1], L2[1]))
    for i = 2:N
        expr = :(_add($expr, _mul(L1[$i], L2[$i])))
    end
    return expr
end

# ──────────────────────────────────────────────
# Cascade: ▷
# ──────────────────────────────────────────────

"""
    ▷(G1::SLH{N}, G2::SLH{N}) where N

Cascade two SLH triples:

``G_1 \\triangleright G_2 = (S_2 S_1,\\; L_2 + S_2 L_1,\\; H_1 + H_2 - \\tfrac{i}{2}(L_2^\\dagger S_2 L_1 - L_1^\\dagger S_2^\\dagger L_2))``

Unicode `\\triangleright<tab>`. See also [`cascade`](@ref).
"""
function ▷(G1::SLH{N}, G2::SLH{N}) where {N}
    S1, L1, H1 = scattering(G1), lindblad(G1), hamiltonian(G1)
    S2, L2, H2 = scattering(G2), lindblad(G2), hamiltonian(G2)

    S_t = _post.(S2 * S1)
    S2L1 = _slh_matvec(S2, L1)
    L_t = SVector{N}(ntuple(i -> _post(_add(L2[i], S2L1[i])), Val(N)))

    L2_adj = SVector{N}(ntuple(i -> _adj(L2[i]), Val(N)))
    cross1 = _slh_dot(L2_adj, S2L1)
    X = _mul(-1im / 2, cross1)

    H_t = _post(_add(_add(H1, H2), _add(X, _adj(X))))

    op_hint = _op_type(G1)
    op_hint === nothing && (op_hint = _op_type(G2))
    return _build_slh(S_t, L_t, H_t, op_hint)
end

▷(a::SLH, b::SLH, c::SLH...) = ▷(a ▷ b, c...)

"""
    cascade(G::SLH...)

Cascade SLH triples from first to last. Alias for [`▷`](@ref).
"""
cascade(args...) = ▷(args...)

# ──────────────────────────────────────────────
# Concatenate: ⊞
# ──────────────────────────────────────────────

"""
    ⊞(G1::SLH{N1}, G2::SLH{N2})

Concatenate (parallel composition) of two SLH triples:

``G_1 \\boxplus G_2 = \\left(\\begin{pmatrix} S_1 & 0 \\\\ 0 & S_2 \\end{pmatrix},\\;
\\begin{pmatrix} L_1 \\\\ L_2 \\end{pmatrix},\\; H_1 + H_2\\right)``

Unicode `\\boxplus<tab>`. See also [`concatenate`](@ref).
"""
@generated function ⊞(G1::SLH{N1}, G2::SLH{N2}) where {N1,N2}
    N = N1 + N2
    s_exprs = []
    for j = 1:N, i = 1:N  # column-major for SMatrix
        if i <= N1 && j <= N1
            push!(s_exprs, :(S1[$i, $j]))
        elseif i > N1 && j > N1
            push!(s_exprs, :(S2[$(i-N1), $(j-N1)]))
        else
            push!(s_exprs, 0)
        end
    end

    quote
        S1, L1, H1 = scattering(G1), lindblad(G1), hamiltonian(G1)
        S2, L2, H2 = scattering(G2), lindblad(G2), hamiltonian(G2)
        S_t = SMatrix{$N,$N}($(s_exprs...))
        L_t = vcat(L1, L2)
        H_t = _add(H1, H2)
        op_hint = _op_type(G1)
        op_hint === nothing && (op_hint = _op_type(G2))
        return _build_slh(S_t, L_t, H_t, op_hint)
    end
end

⊞(a::SLH, b::SLH, c::SLH...) = ⊞(a ⊞ b, c...)

"""
    concatenate(G::SLH...)

Concatenate (parallel composition) of SLH triples. Alias for [`⊞`](@ref).
"""
concatenate(args...) = ⊞(args...)

# ──────────────────────────────────────────────
# Feedback reduction
# dispatch through Val(N-1) so M is a type parameter
# all ntuple calls use Val
# ──────────────────────────────────────────────

function _drop_row_col(S::SMatrix{N,N}, row::Int, col::Int, ::Val{M}) where {N,M}
    SMatrix{M,M}(ntuple(Val(M * M)) do k
        i, j = divrem(k - 1, M) .+ (1, 1)
        ri = i >= row ? i + 1 : i
        cj = j >= col ? j + 1 : j
        S[ri, cj]
    end)
end

function _drop_index(L::SVector{N}, idx::Int, ::Val{M}) where {N,M}
    SVector{M}(ntuple(i -> L[i >= idx ? i + 1 : i], Val(M)))
end

function _get_col_dropped_row(
    S::SMatrix{N,N},
    col::Int,
    drop_row::Int,
    ::Val{M},
) where {N,M}
    SVector{M}(ntuple(i -> S[i >= drop_row ? i + 1 : i, col], Val(M)))
end

function _get_row_dropped_col(
    S::SMatrix{N,N},
    row::Int,
    drop_col::Int,
    ::Val{M},
) where {N,M}
    SVector{M}(ntuple(j -> S[row, j >= drop_col ? j + 1 : j], Val(M)))
end

"""
    feedback(G::SLH{N}, x::Int, y::Int) where N

Apply the SLH feedback reduction rule: connect output port `x` to input port `y`.
Returns `SLH{N-1}`.

See also [`SLH`](@ref), [`▷`](@ref), [`⊞`](@ref).
"""
function feedback(G::SLH{N}, x::Int, y::Int) where {N}
    _feedback_impl(G, x, y, Val(N - 1))
end

function _feedback_impl(G::SLH{N}, x::Int, y::Int, ::Val{M}) where {N,M}
    S = scattering(G)
    L = lindblad(G)
    H = hamiltonian(G)

    @assert 1 <= x <= N && 1 <= y <= N

    valM = Val(M)
    S_xy = S[x, y]
    loop_gain = 1 / (1 - S_xy)

    S_bar = _drop_row_col(S, x, y, valM)
    S_col_y_no_x = _get_col_dropped_row(S, y, x, valM)
    S_row_x_no_y = _get_row_dropped_col(S, x, y, valM)

    S_update = SMatrix{M,M}(ntuple(Val(M * M)) do k
        i, j = divrem(k - 1, M) .+ (1, 1)
        _post(_mul(_mul(S_col_y_no_x[i], loop_gain), S_row_x_no_y[j]))
    end)
    S_red = SMatrix{M,M}(ntuple(Val(M * M)) do k
        i, j = divrem(k - 1, M) .+ (1, 1)
        _post(_add(S_bar[i, j], S_update[i, j]))
    end)

    L_bar = _drop_index(L, x, valM)
    L_x = L[x]
    L_update =
        SVector{M}(ntuple(i -> _post(_mul(_mul(S_col_y_no_x[i], loop_gain), L_x)), valM))
    L_red = SVector{M}(ntuple(i -> _post(_add(L_bar[i], L_update[i])), valM))

    S_col_y_full = SVector{N}(ntuple(i -> S[i, y], Val(N)))
    L_adj = SVector{N}(ntuple(i -> _adj(L[i]), Val(N)))
    term = _mul(_slh_dot(L_adj, S_col_y_full), _mul(loop_gain, L_x))
    H_red = _post(_add(H, _mul(1 / (2im), _add(term, _mul(-1, _adj(term))))))

    return _build_slh(S_red, L_red, H_red, _op_type(G))
end

feedback(G::SLH, connection::Pair{Int,Int}) =
    feedback(G, first(connection), last(connection))

function feedback(G::SLH, connections::Pair{Int,Int}...)
    n = size(scattering(G), 1)
    G_red = G
    for (x, y) in _feedback_maps(n, connections)
        G_red = feedback(G_red, x, y)
    end
    return G_red
end

function _feedback_maps(n::Int, connections)
    output_ports = collect(1:n)
    input_ports = collect(1:n)
    mapped = Tuple{Int,Int}[]

    for connection in connections
        x = first(connection)
        y = last(connection)
        x_now = findfirst(==(x), output_ports)
        y_now = findfirst(==(y), input_ports)
        @assert x_now !== nothing "output port $x has already been eliminated"
        @assert y_now !== nothing "input port $y has already been eliminated"
        push!(mapped, (x_now, y_now))
        deleteat!(output_ports, x_now)
        deleteat!(input_ports, y_now)
    end
    return mapped
end

# ──────────────────────────────────────────────
# Numeric translation of an SLH object
# ──────────────────────────────────────────────

"""
    to_numeric(G::SLH, b::QuantumOpticsBase.Basis; kwargs...)

Translate the Hamiltonian and Lindblad operators of an SLH object `G` into numeric
[QuantumOptics.jl](https://github.com/qojulia/QuantumOptics.jl) operators on the basis `b`.
Returns the tuple `(H_QO, L_QO)`, where `L_QO` is a vector holding one translated operator
per jump operator in `lindblad(G)`. All keyword arguments (`parameter`, `time_parameter`,
`operators`, `adjoint_ops`, `op_type`) are forwarded to
[`SecondQuantizedAlgebra.to_numeric`](@ref).
"""
function SQA.to_numeric(G::SLH, b::QuantumOpticsBase.Basis; kwargs...)
    H_QO = SQA.to_numeric(hamiltonian(G), b; kwargs...)
    L_QO = [SQA.to_numeric(L_, b; kwargs...) for L_ in lindblad(G)]
    return H_QO, L_QO
end
