"""
    quantum_fisher_information(ρ, dρ; regularization=1e-12)

Compute the quantum Fisher information from a density matrix `ρ` and its
parameter derivative `dρ`.

The implementation uses the symmetric logarithmic derivative `L`, defined by
`dρ = (ρL + Lρ) / 2`, and returns `real(tr(dρ * L))`. A small diagonal
`regularization` is added to `ρ` before solving the Sylvester equation so that
rank-deficient states are numerically well behaved.
"""
function quantum_fisher_information(
    ρ::AbstractMatrix,
    dρ::AbstractMatrix;
    regularization = 1e-12,
)
    size(ρ, 1) == size(ρ, 2) || throw(ArgumentError("`ρ` must be square."))
    size(dρ) == size(ρ) || throw(ArgumentError("`dρ` must have the same size as `ρ`."))

    ρ_reg = ρ + regularization * LinearAlgebra.I
    L = LinearAlgebra.sylvester(ρ_reg, ρ_reg, -2 * dρ)
    return real(LinearAlgebra.tr(dρ * L))
end

function quantum_fisher_information(
    ρ::QuantumOpticsBase.AbstractOperator,
    dρ::QuantumOpticsBase.AbstractOperator;
    kwargs...,
)
    return quantum_fisher_information(Matrix(ρ.data), Matrix(dρ.data); kwargs...)
end

_matrix_data(A::AbstractMatrix) = Matrix(A)
_matrix_data(A::QuantumOpticsBase.AbstractOperator) = Matrix(A.data)

function _assert_square_same_size(A, B; names = ("first argument", "second argument"))
    size(A, 1) == size(A, 2) || throw(ArgumentError("`$(names[1])` must be square."))
    size(B) == size(A) || throw(ArgumentError("`$(names[2])` must have the same size as `$(names[1])`."))
    return nothing
end

"""
    povm_probabilities(ρ, measurement)

Return the probabilities ``p_i = Tr[M_i ρ]`` for a POVM or projective
measurement `measurement = [M₁, M₂, ...]`.

`ρ` and each measurement element may be matrices or `QuantumOpticsBase`
operators.
"""
function povm_probabilities(ρ, measurement)
    ρm = _matrix_data(ρ)
    return [real(LinearAlgebra.tr(_matrix_data(M) * ρm)) for M in measurement]
end

"""
    classical_fisher_information(ρ, dρ, measurement; probability_floor=1e-12,
                                 derivative_floor=1e-12)

Compute the classical Fisher information for a concrete POVM or projective
measurement `measurement`.

For ``p_i = Tr[M_i ρ]`` and ``∂p_i = Tr[M_i ∂ρ]``, this returns

```math
F_C = \\sum_i \\frac{(\\partial_\\theta p_i)^2}{p_i}.
```

Outcomes with probability below `probability_floor` and derivative below
`derivative_floor` are ignored. If an impossible outcome has a nonzero
derivative, the Fisher information is infinite.
"""
function classical_fisher_information(
    ρ,
    dρ,
    measurement;
    probability_floor = 1e-12,
    derivative_floor = 1e-12,
)
    ρm = _matrix_data(ρ)
    dρm = _matrix_data(dρ)
    _assert_square_same_size(ρm, dρm; names = ("ρ", "dρ"))

    F = 0.0
    for M in measurement
        Mm = _matrix_data(M)
        _assert_square_same_size(ρm, Mm; names = ("ρ", "measurement element"))
        p = real(LinearAlgebra.tr(Mm * ρm))
        dp = real(LinearAlgebra.tr(Mm * dρm))

        if p < -probability_floor
            throw(DomainError(p, "measurement probabilities must be nonnegative"))
        elseif p <= probability_floor
            abs(dp) <= derivative_floor && continue
            return Inf
        end
        F += abs2(dp) / p
    end
    return F
end

"""
    projective_measurement(A; atol=1e-10)

Return the projectors onto the eigenspaces of a Hermitian observable `A`.
Degenerate eigenvalues within `atol` are grouped into one projector.
"""
function projective_measurement(A::AbstractMatrix; atol = 1e-10)
    isapprox(A, A'; atol) || throw(ArgumentError("observable must be Hermitian."))
    F = LinearAlgebra.eigen(LinearAlgebra.Hermitian(Matrix(A)))
    projectors = Matrix{eltype(F.vectors)}[]
    used = falses(length(F.values))

    for i in eachindex(F.values)
        used[i] && continue
        inds = findall(j -> !used[j] && abs(F.values[j] - F.values[i]) <= atol, eachindex(F.values))
        P = zeros(eltype(F.vectors), size(A, 1), size(A, 2))
        for j in inds
            v = F.vectors[:, j]
            P += v * v'
            used[j] = true
        end
        push!(projectors, P)
    end
    return projectors
end

function projective_measurement(A::QuantumOpticsBase.AbstractOperator; kwargs...)
    return [
        QuantumOpticsBase.Operator(A.basis_l, A.basis_r, P)
        for P in projective_measurement(Matrix(A.data); kwargs...)
    ]
end

function _split_complex_matrix(v, d)
    half = length(v) ÷ 2
    return reshape(v[1:half] + im * v[(half + 1):end], d, d)
end

function _density_state_vector(ρ)
    ρvec = vec(Matrix(ρ.data))
    return vcat(real.(ρvec), imag.(ρvec))
end

function _evolve_density(H, J, ρ0, t; kwargs...)
    iszero(t) && return ρ0
    _, ρt = timeevolution.master((zero(float(real(t))), t), ρ0, H, J; saveat = [t], save_everystep = false, kwargs...)
    return last(ρt)
end

"""
    parameter_derivative(G::SLH, b, ρ0, t; estimate, parameter=Dict(), operators=Dict(), kwargs...)

Evolve the SLH model `G` to time `t` and return `(ρ, dρ)`, where `dρ` is the
derivative of the density matrix with respect to the scalar symbolic parameter
`estimate`.

The keyword `parameter` gives the numerical working point and must include
`estimate`. The Hamiltonian and Lindblad operators are obtained with
[`translate_qo`](@ref), so `estimate` may appear in `H` or `L`. The derivative is
computed with `ForwardDiff` through the `QuantumOptics.timeevolution.master`
solve. Extra keyword arguments are forwarded to the time-evolution solver.
"""
function parameter_derivative(
    G::SLH,
    b::QuantumOpticsBase.Basis,
    ρ0::QuantumOpticsBase.AbstractOperator,
    t;
    estimate,
    parameter = Dict(),
    operators = Dict(),
    kwargs...,
)
    haskey(parameter, estimate) || throw(ArgumentError("`parameter` must contain the estimated parameter `$estimate`."))
    θ0 = parameter[estimate]
    θ0 isa Real || throw(ArgumentError("`estimate` must have a real scalar value, got `$θ0`."))
    d = size(ρ0.data, 1)

    function state_vector(θ)
        p = Dict{Any,Any}(parameter)
        p[estimate] = θ[1]
        H, J = translate_qo(G, b; parameter = p, operators)
        return _density_state_vector(_evolve_density(H, J, ρ0, t; kwargs...))
    end

    θ = [θ0]
    ρvec = state_vector(θ)
    dρvec = vec(ForwardDiff.jacobian(state_vector, θ))
    ρ = _split_complex_matrix(ρvec, d)
    dρ = _split_complex_matrix(dρvec, d)
    return (
        QuantumOpticsBase.Operator(ρ0.basis_l, ρ0.basis_r, ρ),
        QuantumOpticsBase.Operator(ρ0.basis_l, ρ0.basis_r, dρ),
    )
end
