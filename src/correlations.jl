"""
    correlation_matrix(T, ρt, f, Ls; kwargs...)
    correlation_matrix(T, ρt, H, J, Ls; kwargs...)

Compute the two-time correlation matrix
``g^{(1)}(t_1, t_2) = \\langle L_s^\\dagger(t_1) L_s(t_2) \\rangle``
on the time grid `T`. Writes directly into output matrix.
"""
function correlation_matrix(T::Vector, ρt::Vector, f::Function, Ls::Function; kwargs...)
    Ls_vec = Ls.(T)
    Ls_dag_vec = dagger.(Ls_vec)
    _correlation_loop(T, ρt, Ls_vec, Ls_dag_vec) do T_slice, ρ0
        timeevolution.master_dynamic(T_slice, ρ0, f; kwargs...)
    end
end

function correlation_matrix(T::Vector, ρt::Vector, f::Function, Ls; kwargs...)
    Ls_dag = dagger(Ls)
    l_T = length(T)
    Ls_vec = fill(Ls, l_T)
    Ls_dag_vec = fill(Ls_dag, l_T)
    _correlation_loop(T, ρt, Ls_vec, Ls_dag_vec) do T_slice, ρ0
        timeevolution.master_dynamic(T_slice, ρ0, f; kwargs...)
    end
end

function correlation_matrix(T::Vector, ρt::Vector, H, J::Vector, Ls; kwargs...)
    Ls_dag = dagger(Ls)
    l_T = length(T)
    Ls_vec = fill(Ls, l_T)
    Ls_dag_vec = fill(Ls_dag, l_T)
    _correlation_loop(T, ρt, Ls_vec, Ls_dag_vec) do T_slice, ρ0
        timeevolution.master(T_slice, ρ0, H, J; kwargs...)
    end
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
