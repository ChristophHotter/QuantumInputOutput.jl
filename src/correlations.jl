"""
    correlation_matrix(T, ρt, f, Js; kwargs...)
    correlation_matrix(T, ρt, H, J, Js; kwargs...)

Compute the two-time correlation matrix
``g^{(1)}(t_1, t_2) = \\langle J_s^\\dagger(t_1) J_s(t_2) \\rangle``
on the time grid `T`. Writes directly into output matrix.

Supply the dynamics either as a `master_dynamic`-style function `f(t, ρ)`, or as operators
passed straight to the solver: a time-dependent `H` (e.g. the `TimeDependentSum` from
[`to_numeric`](@ref)) with jump operators `J`, or a constant `H` with constant `J`. The
operator form is much faster for time-dependent problems (the integrator is built once).
`Js` is either a constant operator or a function `Js(t)` returning the operator at `t`.

The returned matrix is a `Hermitian` wrapper. To extract the dominant temporal modes,
diagonalize it with `eigen(g1_m)`. When only the leading modes are needed, te cheaper eigenvalue-range method can be used, e.g. `eigen(g1_m, (n-4):n)` for the five dominant modes, where `n = size(g1_m, 1)`.
"""
function correlation_matrix(T::Vector, ρt::Vector, f::Function, Js; kwargs...)
    Js_vec, Js_dag_vec = _sample_operator_and_adjoint(T, Js)
    _correlation_loop(T, ρt, Js_vec, Js_dag_vec) do T_slice, ρ0
        timeevolution.master_dynamic(T_slice, ρ0, f; kwargs...)
    end
end

function correlation_matrix(
    T::Vector,
    ρt::Vector,
    H::QuantumOpticsBase.AbstractTimeDependentOperator,
    J::AbstractVector,
    Js;
    kwargs...,
)
    Js_vec, Js_dag_vec = _sample_operator_and_adjoint(T, Js)
    _correlation_loop(T, ρt, Js_vec, Js_dag_vec) do T_slice, ρ0
        timeevolution.master_dynamic(T_slice, ρ0, copy(H), [copy(j) for j in J]; kwargs...)
    end
end

function correlation_matrix(T::Vector, ρt::Vector, H, J::AbstractVector, Js; kwargs...)
    Js_vec, Js_dag_vec = _sample_operator_and_adjoint(T, Js)
    _correlation_loop(T, ρt, Js_vec, Js_dag_vec) do T_slice, ρ0
        timeevolution.master(T_slice, ρ0, H, J; kwargs...)
    end
end

function _sample_operator_and_adjoint(T::Vector, op::Function)
    vals = op.(T)
    return vals, dagger.(vals)
end
function _sample_operator_and_adjoint(
    ::Vector,
    op::QuantumOpticsBase.AbstractTimeDependentOperator,
)
    throw(
        ArgumentError(
            "`Js` is a time-dependent operator ($(nameof(typeof(op)))); pass it as a " *
            "function `Js(t)` returning the concrete operator at time `t`.",
        ),
    )
end
function _sample_operator_and_adjoint(T::Vector, op)
    nt = length(T)
    return fill(op, nt), fill(dagger(op), nt)
end

function _correlation_loop(solve_fn, T, ρt, Js_vec, Js_dag_vec)
    l_T = length(T)
    @assert l_T == length(ρt)

    g1_m = zeros(ComplexF64, l_T, l_T)
    # Each iteration solves an independent master equation — parallelise
    Threads.@threads for it = 1:(l_T-1)
        ρ0_it = Js_vec[it] * ρt[it]
        τ_, ρ_bar_τ = solve_fn(@view(T[it:end]), ρ0_it)

        @inbounds for i in eachindex(ρ_bar_τ)
            val = expect(Js_dag_vec[it+i-1], ρ_bar_τ[i])
            g1_m[it, it+i-1] = val
            g1_m[it+i-1, it] = conj(val)
        end
    end
    g1_m[l_T, l_T] = expect(Js_dag_vec[l_T], Js_vec[l_T] * ρt[l_T])
    return Hermitian(g1_m)
end
