"""
    two_time_corr_matrix(T, ρt, f, Ls; kwargs...)
    two_time_corr_matrix(T, ρt, H, J, Ls; kwargs...)

Compute the two-time correlation matrix ``g^{(1)}(t_1, t_2) = \\langle L_s^\\dagger(t_1) L_s(t_2) \\rangle`` 
on the time grid `T` for the operator `Ls`.
The first method supports time-dependent generators with the function `f`; the second is for time-independent `H` and `J`.
"""
function two_time_corr_matrix(T::Vector, ρt::Vector, f::Function, Ls::Function; kwargs...)
    l_T = length(T)
    @assert l_T == length(ρt)
    Ls_ls = Ls.(T)
    Ls_ls_dag = dagger.(Ls_ls)
    ρ0_ = [Ls_ls[i] * ρt[i] for i = 1:l_T]

    g1_m = zeros(ComplexF64, l_T, l_T)
    for it = 1:(l_T-1)
        τ_, ρ_bar_τ = timeevolution.master_dynamic(T[it:end], ρ0_[it], f; kwargs...)

        g1 = [expect(Ls_ls_dag[it+i-1], ρ_bar_τ[i]) for i = 1:length(τ_)]
        g1_m[it, it:end] = g1
        g1_m[it:end, it] = adjoint.(g1)
    end
    return g1_m
end

# time-dependent problem but constant Ls
function two_time_corr_matrix(T::Vector, ρt::Vector, f::Function, Ls; kwargs...)
    l_T = length(T)
    @assert l_T == length(ρt)
    Ls_dag = dagger.(Ls)
    ρ0_ = [Ls * ρt[i] for i = 1:l_T]

    g1_m = zeros(ComplexF64, l_T, l_T)
    for it = 1:(l_T-1)
        τ_, ρ_bar_τ = timeevolution.master_dynamic(T[it:end], ρ0_[it], f; kwargs...)

        g1 = [expect(Ls_dag, ρ_bar_τ[i]) for i = 1:length(τ_)]
        g1_m[it, it:end] = g1
        g1_m[it:end, it] = adjoint.(g1)
    end
    return g1_m
end

# two_time_corr_matrix for time-independent problems
function two_time_corr_matrix(T::Vector, ρt::Vector, H, J::Vector, Ls; kwargs...)
    l_T = length(T)
    @assert l_T == length(ρt)
    Ls_dag = dagger(Ls)
    ρ0_ = [Ls * ρt[i] for i = 1:l_T]

    g1_m = zeros(ComplexF64, l_T, l_T)
    for it = 1:(l_T-1)
        τ_, ρ_bar_τ = timeevolution.master(T[it:end], ρ0_[it], H, J; kwargs...)

        g1 = [expect(Ls_dag, ρ_bar_τ[i]) for i = 1:length(τ_)]
        g1_m[it, it:end] = g1
        g1_m[it:end, it] = adjoint.(g1)
    end
    return g1_m
end
