"""
    correlation_matrix(T, ρt, f, Ls; kwargs...)
    correlation_matrix(T, ρt, H, J, Ls; kwargs...)

Compute the two-time correlation matrix
``g^{(1)}(t_1, t_2) = \\langle L_s^\\dagger(t_1) L_s(t_2) \\rangle``
on the time grid `T`. Writes directly into output matrix.
"""
function correlation_matrix(T::Vector, ρt::Vector, f::Function, Ls::Function; kwargs...)
    l_T = length(T)
    @assert l_T == length(ρt)
    Ls_ls = Ls.(T)
    Ls_ls_dag = dagger.(Ls_ls)

    g1_m = zeros(ComplexF64, l_T, l_T)
    for it = 1:(l_T-1)
        ρ0_it = Ls_ls[it] * ρt[it]
        τ_, ρ_bar_τ = timeevolution.master_dynamic(T[it:end], ρ0_it, f; kwargs...)

        @inbounds for i in eachindex(ρ_bar_τ)
            val = expect(Ls_ls_dag[it+i-1], ρ_bar_τ[i])
            g1_m[it, it+i-1] = val
            g1_m[it+i-1, it] = conj(val)
        end
    end
    return g1_m
end

function correlation_matrix(T::Vector, ρt::Vector, f::Function, Ls; kwargs...)
    l_T = length(T)
    @assert l_T == length(ρt)
    Ls_dag = dagger(Ls)

    g1_m = zeros(ComplexF64, l_T, l_T)
    for it = 1:(l_T-1)
        ρ0_it = Ls * ρt[it]
        τ_, ρ_bar_τ = timeevolution.master_dynamic(T[it:end], ρ0_it, f; kwargs...)

        @inbounds for i in eachindex(ρ_bar_τ)
            val = expect(Ls_dag, ρ_bar_τ[i])
            g1_m[it, it+i-1] = val
            g1_m[it+i-1, it] = conj(val)
        end
    end
    return g1_m
end

function correlation_matrix(T::Vector, ρt::Vector, H, J::Vector, Ls; kwargs...)
    l_T = length(T)
    @assert l_T == length(ρt)
    Ls_dag = dagger(Ls)

    g1_m = zeros(ComplexF64, l_T, l_T)
    for it = 1:(l_T-1)
        ρ0_it = Ls * ρt[it]
        τ_, ρ_bar_τ = timeevolution.master(T[it:end], ρ0_it, H, J; kwargs...)

        @inbounds for i in eachindex(ρ_bar_τ)
            val = expect(Ls_dag, ρ_bar_τ[i])
            g1_m[it, it+i-1] = val
            g1_m[it+i-1, it] = conj(val)
        end
    end
    return g1_m
end
