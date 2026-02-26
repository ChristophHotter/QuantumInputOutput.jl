# creates matrix of the two time correlation function
function two_time_corr_matrix(T_ls::Vector, ρt::Vector, f::Function, Ls::Function; abstol=1e-6, reltol=1e-6) #kwargs...) TODO
    l_T_ls = length(T_ls)
    @assert l_T_ls == length(ρt)
    Ls_ls = Ls.(T_ls)
    Ls_ls_dag = dagger.(Ls_ls)
    ρ0_ = [Ls_ls[i] * ρt[i] for i = 1:l_T_ls]

    g1_m = zeros(ComplexF64, l_T_ls, l_T_ls)
    for it = 1:l_T_ls-1
        τ_, ρ_bar_τ = timeevolution.master_dynamic(T_ls[it:end], ρ0_[it], f; abstol = abstol, reltol = reltol) #kwargs...) TODO

        g1 = [expect(Ls_ls_dag[it + i - 1], ρ_bar_τ[i]) for i = 1:length(τ_)]
        g1_m[it, it:end] = g1
        g1_m[it:end, it] = adjoint.(g1)
    end
    return g1_m
end

# two_time_corr_matrix for time-independent problems
function two_time_corr_matrix(T_ls::Vector, ρt::Vector, H, J::Vector, Ls; kwargs...)
    l_T_ls = length(T_ls)
    @assert l_T_ls == length(ρt)
    Ls_dag = dagger(Ls)
    ρ0_ = [Ls * ρt[i] for i = 1:l_T_ls]

    g1_m = zeros(ComplexF64, l_T_ls, l_T_ls)
    for it = 1:l_T_ls-1
        τ_, ρ_bar_τ = timeevolution.master(T_ls[it:end], ρ0_[it], H, J; kwargs...)

        g1 = [expect(Ls_dag, ρ_bar_τ[i]) for i = 1:length(τ_)]
        g1_m[it, it:end] = g1
        g1_m[it:end, it] = adjoint.(g1)
    end
    return g1_m
end
