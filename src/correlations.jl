"""
    correlation_matrix(T, ρt, f, Ls; kwargs...)
    correlation_matrix(T, ρt, H, J, Ls; kwargs...)

Compute the two-time correlation matrix
``g^{(1)}(t_1, t_2) = \\langle L_s^\\dagger(t_1) L_s(t_2) \\rangle``
on the time grid `T`. Writes directly into output matrix.

The dynamics can be supplied either as a `master_dynamic`-style function `f(t, ρ)` or,
much more efficiently for time-dependent problems, as operators passed straight to the
solver: a time-dependent `H` (an `AbstractTimeDependentOperator`, e.g. the
`TimeDependentSum` returned by [`to_numeric`](@ref)) together with its jump operators `J`,
or a constant `H` with constant `J`. Passing the operators directly lets the solver build
the integrator once instead of rebuilding it from `f`'s return value at every step.

`Ls` may be a constant operator or a function `Ls(t)` returning the (concrete) operator at
time `t`.
"""
function correlation_matrix(T::Vector, ρt::Vector, f::Function, Ls; kwargs...)
    Ls_vec, Ls_dag_vec = _ls_vectors(T, Ls)
    _correlation_loop(T, ρt, Ls_vec, Ls_dag_vec) do T_slice, ρ0
        timeevolution.master_dynamic(T_slice, ρ0, f; kwargs...)
    end
end

function correlation_matrix(
    T::Vector,
    ρt::Vector,
    H::QuantumOpticsBase.AbstractTimeDependentOperator,
    J::AbstractVector,
    Ls;
    kwargs...,
)
    Ls_vec, Ls_dag_vec = _ls_vectors(T, Ls)
    # `master_dynamic` mutates the operator's current time via `set_time!`, so each parallel
    # iteration solves with its own copy to avoid a data race on the shared `H`/`J`.
    _correlation_loop(T, ρt, Ls_vec, Ls_dag_vec) do T_slice, ρ0
        timeevolution.master_dynamic(T_slice, ρ0, copy(H), [copy(j) for j in J]; kwargs...)
    end
end

function correlation_matrix(T::Vector, ρt::Vector, H, J::AbstractVector, Ls; kwargs...)
    Ls_vec, Ls_dag_vec = _ls_vectors(T, Ls)
    _correlation_loop(T, ρt, Ls_vec, Ls_dag_vec) do T_slice, ρ0
        timeevolution.master(T_slice, ρ0, H, J; kwargs...)
    end
end

# Build the per-time `Ls`/`Ls†` vectors. A function `Ls(t)` is sampled on `T`; a constant
# operator is broadcast to every time.
function _ls_vectors(T::Vector, Ls::Function)
    Ls_vec = Ls.(T)
    return Ls_vec, dagger.(Ls_vec)
end
function _ls_vectors(T::Vector, Ls)
    l_T = length(T)
    return fill(Ls, l_T), fill(dagger(Ls), l_T)
end

function _correlation_loop(solve_fn, T, ρt, Ls_vec, Ls_dag_vec)
    l_T = length(T)
    @assert l_T == length(ρt)

    g1_m = zeros(ComplexF64, l_T, l_T)
    # Each iteration solves an independent master equation — parallelise
    Threads.@threads for it = 1:(l_T-1)
        ρ0_it = Ls_vec[it] * ρt[it]
        τ_, ρ_bar_τ = solve_fn(@view(T[it:end]), ρ0_it)

        @inbounds for i in eachindex(ρ_bar_τ)
            val = expect(Ls_dag_vec[it+i-1], ρ_bar_τ[i])
            g1_m[it, it+i-1] = val
            g1_m[it+i-1, it] = conj(val)
        end
    end
    return g1_m
end
