"""
    SLH

SLH triple with scattering matrix `S`, Lindblad term `L` and Hamiltonian `H`. 
`S` and `L` can also be vectors of scattering matrices and Linblad terms. 

See also [`▷`](@ref) and [`⊞`](@ref)
"""
struct SLH{T,LT,H,S<:AbstractMatrix{T},L<:AbstractVector{LT}}
    scattering::S
    lindblad::L
    hamiltonian::H

    function SLH(
        scattering::S,
        lindblad::L,
        hamiltonian::H,
    ) where {T,LT,H,S<:AbstractMatrix{T},L<:AbstractVector{LT}}

        @assert size(scattering, 1) == length(lindblad)
        new{T,LT,H,S,L}(scattering, lindblad, hamiltonian)
    end
end

# -----------------------
# QuantumOptics operators
# -----------------------

"""
    SLHqo

SLH triple where `L` and `H` are QuantumOptics.jl operators or time-dependent functions
returning such operators. This is useful when you build models directly in a numeric basis,
without symbolic translation.
"""
struct SLHqo{S<:AbstractMatrix,L<:AbstractVector,H}
    scattering::S
    lindblad::L
    hamiltonian::H

    function SLHqo(
        scattering::S,
        lindblad::L,
        hamiltonian::H,
    ) where {S<:AbstractMatrix,L<:AbstractVector,H}
        @assert size(scattering, 1) == length(lindblad)
        if !(hamiltonian isa Function || hamiltonian isa QuantumOpticsBase.AbstractOperator)
            error(
                "SLHqo expects H to be a QuantumOptics operator or a time-dependent function returning one; got $(typeof(hamiltonian))",
            )
        end
        for (i, l) in pairs(lindblad)
            if !(l isa Function || l isa QuantumOpticsBase.AbstractOperator)
                error(
                    "SLHqo expects L[$i] to be a QuantumOptics operator or a time-dependent function returning one; got $(typeof(l))",
                )
            end
        end
        new{S,L,H}(scattering, lindblad, hamiltonian)
    end
end

function SLHqo(S, L::AbstractVector, H)
    S_ = Matrix(I, length(L), length(L)) * S
    return SLHqo(S_, L, H)
end
function SLHqo(S, L, H)
    S_ = ones(1, 1) * S
    L_ = [L]
    return SLHqo(S_, L_, H)
end

"""
    get_scattering(G)

Return the scattering matrix `S` of an SLH or SLHqo object.
"""
get_scattering(slh::SLHqo) = slh.scattering
"""
    get_lindblad(G)

Return the Lindblad vector `L` of an SLH or SLHqo object.
"""
get_lindblad(slh::SLHqo) = slh.lindblad
"""
    get_hamiltonian(G)

Return the Hamiltonian `H` of an SLH or SLHqo object.
"""
get_hamiltonian(slh::SLHqo) = slh.hamiltonian

function Base.isequal(slh1::SLHqo, slh2::SLHqo)
    isequal(get_scattering(slh1), get_scattering(slh2)) &&
        isequal(get_lindblad(slh1), get_lindblad(slh2)) &&
        isequal(get_hamiltonian(slh1), get_hamiltonian(slh2))
end

_is_time_dep(x) = x isa Function
_to_func(x) = x isa Function ? x : (t -> x)

# SLH(scattering::S, lindblad::L, hamiltonian::H) where {S,H,L} = SLH{S,H,L}(scattering,lindblad,hamiltonian)
function SLH(S, L::AbstractVector, H)
    S_ = Matrix{typeof(S)}(I, length(L), length(L))*S
    return SLH(S_, L, H)
end
function SLH(S, L, H)
    S_ = ones(1, 1)*S
    L_ = [L]
    return SLH(S_, L_, H)
end

function Base.isequal(slh1::SLH, slh2::SLH)
    isequal(get_scattering(slh1), get_scattering(slh2)) &&
        isequal(get_lindblad(slh1), get_lindblad(slh2)) &&
        isequal(get_hamiltonian(slh1), get_hamiltonian(slh2))
end

get_scattering(slh::SLH) = slh.scattering
get_lindblad(slh::SLH) = slh.lindblad
get_hamiltonian(slh::SLH) = slh.hamiltonian

"""
    ▷(G::SLH...)

Creates a new SLH triple by cascading the SLH triples from first to last according to 
the rule 

``SLH_1 \\triangleright SLH_2 = ( S_2 S1, L_2 + S_2 L_1, H_1 + H_2 - \\frac{i}{2} L_2^\\dagger S_2 L_1 - L_1^\\dagger S_2^\\dagger L_2 ) ``

Unicode `\\triangleright<tab>` alias of [`cascade`](@ref)
"""
function ▷(G1::SLH, G2::SLH) #\triangleright
    # function cascade(G1::SLH,G2::SLH) # G1 ▷ G2 = G2 ◁ G1
    S1 = get_scattering(G1);
    L1 = get_lindblad(G1);
    H1 = get_hamiltonian(G1)
    S2 = get_scattering(G2);
    L2 = get_lindblad(G2);
    H2 = get_hamiltonian(G2)

    lL1 = length(L1)
    lL2 = length(L2)
    @assert lL1==lL2

    S_t = simplify.(S2*S1)
    L_t = simplify.(L2 + S2*L1)
    # H_t = (H1 + H2 - 1im/2*(L2'S2*L1 - L1'S2'L2)) # the adjoint also creates the transpose! 
    # H_t = (H1 + H2 - 1im/2*(SQA._adjoint(L2)*S2*L1 - SQA._adjoint(L1)*SQA._adjoint(S2)*L2))
    L1_ct = SQA._adjoint.(L1)
    L2_ct = SQA._adjoint.(L2)
    S2_ct = SQA._adjoint.(S2)

    H_t = simplify(
        H1 + H2 -
        1im/2*(
            sum(L2_ct[it]*(S2*L1)[it] for it = 1:lL1) -
            sum(L1_ct[it]*(S2_ct*L2)[it] for it = 1:lL1)
        ),
    )
    return SLH(S_t, L_t, H_t)
end
▷(a::SLH, b::SLH, c::SLH...) = ▷(a▷b, c...)
# G's in the order as they are cascaded (from left to right)

"""
    ▷(G::SLHqo...)

Cascades SLHqo triples from first to last, allowing time-dependent `L` or `H`.
If any input `L` or `H` is time-dependent, the result uses time-dependent functions.
"""
function ▷(G1::SLHqo, G2::SLHqo)
    S1 = get_scattering(G1);
    L1 = get_lindblad(G1);
    H1 = get_hamiltonian(G1)
    S2 = get_scattering(G2);
    L2 = get_lindblad(G2);
    H2 = get_hamiltonian(G2)

    lL1 = length(L1)
    lL2 = length(L2)
    @assert lL1 == lL2

    S_t = S2 * S1

    L1f = map(_to_func, L1)
    L2f = map(_to_func, L2)
    H1f = _to_func(H1)
    H2f = _to_func(H2)
    time_dep =
        any(_is_time_dep, L1) ||
        any(_is_time_dep, L2) ||
        _is_time_dep(H1) ||
        _is_time_dep(H2)

    if time_dep
        function L_t(i)
            return t -> begin
                L1_vec = [f(t) for f in L1f]
                L2_vec = [f(t) for f in L2f]
                (L2_vec+S2*L1_vec)[i]
            end
        end
        L_out = [L_t(i) for i = 1:lL1]

        H_out =
            t -> begin
                L1_vec = [f(t) for f in L1f]
                L2_vec = [f(t) for f in L2f]
                term1 = sum(adjoint(L2_vec[i]) * (S2*L1_vec)[i] for i = 1:lL1)
                term2 = sum(adjoint(L1_vec[i]) * (adjoint(S2)*L2_vec)[i] for i = 1:lL1)
                H1f(t) + H2f(t) - 1im / 2 * (term1 - term2)
            end
        return SLHqo(S_t, L_out, H_out)
    end

    L_out = L2 + S2 * L1
    term1 = sum(adjoint(L2[i]) * (S2*L1)[i] for i = 1:lL1)
    term2 = sum(adjoint(L1[i]) * (adjoint(S2)*L2)[i] for i = 1:lL1)
    H_out = H1 + H2 - 1im / 2 * (term1 - term2)
    return SLHqo(S_t, L_out, H_out)
end
▷(a::SLHqo, b::SLHqo, c::SLHqo...) = ▷(a▷b, c...)

"""
    cascade(G::SLH...)

Creates a new SLH triple by cascading the SLH triples from first to last according to 
the rule 

```math
SLH_1 \\triangleright SLH_2 = ( S_2 S1, L_2 + S_2 L_1, H_1 + H_2 - \\frac{i}{2} L_2^\\dagger S_2 L_1 - L_1^\\dagger S_2^\\dagger L_2 ) 
```

See also [`▷`](@ref). 
"""
cascade(args...) = ▷(args...)
# ◁(G1::SLH,G2::SLH) = ▷(G2,G1) # unknown unicode character

"""
    ⊞(G::SLH...)

Creates a new SLH triple by concatenating the SLH triples according to 
the rule 

```math
SLH_1 \\boxplus SLH_2 = \\left( \\begin{pmatrix} S_1 & 0 \\; 0 & S_2 \\end{pmatrix}, \\begin{pmatrix} L_1 \\; L_2 \\end{pmatrix}, H_1 + H_2 \\right)
```

Unicode `\\boxplus<tab>` alias of [`concatenate`](@ref)
"""
function ⊞(G1::SLH, G2::SLH) #\boxplus
    S1 = get_scattering(G1);
    L1 = get_lindblad(G1);
    H1 = get_hamiltonian(G1)
    S2 = get_scattering(G2);
    L2 = get_lindblad(G2);
    H2 = get_hamiltonian(G2)

    lS1 = size(S1, 1);
    lS2 = size(S2, 1)

    S_t = Matrix{Any}(undef, lS1+lS2, lS1+lS2) # that is pretty ugly (but it works)
    S_t .= zeros(lS1+lS2, lS1+lS2)
    S_t[1:lS1, 1:lS1] .= S1
    S_t[(1+lS1):(lS1+lS2), (1+lS1):(lS1+lS2)] .= S2

    L_t = [L1; L2]
    H_t = H1 + H2

    return SLH(S_t, L_t, H_t)
end
⊞(a::SLH, b::SLH, c::SLH...) = ⊞(a⊞b, c...)

"""
    ⊞(G::SLHqo...)

Concatenates SLHqo triples, allowing time-dependent `L` or `H`.
If any input `L` or `H` is time-dependent, the result uses time-dependent functions.
"""
function ⊞(G1::SLHqo, G2::SLHqo)
    S1 = get_scattering(G1);
    L1 = get_lindblad(G1);
    H1 = get_hamiltonian(G1)
    S2 = get_scattering(G2);
    L2 = get_lindblad(G2);
    H2 = get_hamiltonian(G2)

    lS1 = size(S1, 1);
    lS2 = size(S2, 1)
    T_S = promote_type(eltype(S1), eltype(S2))
    S_t = Matrix{T_S}(undef, lS1 + lS2, lS1 + lS2)
    fill!(S_t, 0)
    S_t[1:lS1, 1:lS1] .= S1
    S_t[(1+lS1):(lS1+lS2), (1+lS1):(lS1+lS2)] .= S2

    L1f = map(_to_func, L1)
    L2f = map(_to_func, L2)
    time_dep =
        any(_is_time_dep, L1) ||
        any(_is_time_dep, L2) ||
        _is_time_dep(H1) ||
        _is_time_dep(H2)

    if time_dep
        L_out = vcat(L1f, L2f)
        H_out = t -> _to_func(H1)(t) + _to_func(H2)(t)
        return SLHqo(S_t, L_out, H_out)
    end

    L_out = [L1; L2]
    H_out = H1 + H2
    return SLHqo(S_t, L_out, H_out)
end
⊞(a::SLHqo, b::SLHqo, c::SLHqo...) = ⊞(a⊞b, c...)

"""
    concatenate(G::SLH...)

Creates a new SLH triple by concatenating the SLH triples according to 
the rule 

```math
SLH_1 \\boxplus SLH_2 = \\left( \\begin{pmatrix} S_1 & 0 \\; 0 & S_2 \\end{pmatrix}, \\begin{pmatrix} L_1 \\; L_2 \\end{pmatrix}, H_1 + H_2 \\right)
```

See also [`⊞`](@ref).
"""
concatenate(args...) = ⊞(args...)
concatenation(args...) = concatenate(args...) # alias

Base.length(h::SecondQuantizedAlgebra.ConcreteHilbertSpace) = 1
Base.length(h::ProductSpace) = length(h.spaces)
