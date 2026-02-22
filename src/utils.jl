# move to SQA
function numeric_average(op::QC.QNumber, state::Vector; kwargs...)
    op_num = sparse(to_numeric(op, state[1]; kwargs...))
    # TODO: sparse, dense
    return QuantumOpticsBase.expect(op_num, state)
end
function numeric_average(avg::Average, state::Vector; kwargs...)
    op = undo_average(op)
    return numeric_average(op, state; kwargs...)
end

expect(avg::Average, state; kwargs...) = numeric_average(avg, state; kwargs...)
expect(op::QNumber, state; kwargs...) = numeric_average(op, state; kwargs...)

######################################
### functions for virtual cavities ###
######################################

_tol_div = 1e-10 #12
_extrapolate = ExtrapolationType.Extension
_ϵu = 1e-10
_ϵv = 1e-10

# return g_u(t) from u(t) - single input mode
function u_to_gu(u::Vector,T::Vector)
    l_T = length(T)
    ∫u2_t = cumul_integrate(T, abs2.(u)) .+ 0im
    gu_t = zeros(ComplexF64,l_T)
    for i=2:l_T
        if abs(sqrt(1 - ∫u2_t[i])) > _tol_div
            gu_t[i] = u[i]' / sqrt(1 - ∫u2_t[i]+ _ϵu )
        end #else 0
        # gu_t[i] = u[i]' / max(sqrt(1 - ∫u2_t[i]), _tol_div) # TODO: test
    end
    return LinearInterpolation(gu_t, T; extrapolation=_extrapolate)
end
u_to_gu(u::Function,T::Vector) = u_to_gu(u.(T),T)
u_to_gu(u::LinearInterpolation,T::Vector) = u_to_gu(u.(T),T)

# return g_v(t) from v(t) - single output mode
function v_to_gv(v::Vector,T::Vector)
    l_T = length(T)
    ∫v2_t = cumul_integrate(T, abs2.(v)) .+ 0im
    gv_t = zeros(ComplexF64,l_T)
    for i=2:l_T
        if abs(sqrt( ∫v2_t[i] )) > _tol_div
            gv_t[i] = -v[i]' / sqrt( ∫v2_t[i] + _ϵv) 
        end # else 0
        # gv_t[i] = -v[i]' / max(sqrt( ∫v2_t[i]), _tol_div) # TODO: test
    end
    return LinearInterpolation(gv_t, T; extrapolation=_extrapolate)
end
v_to_gv(v::Function,T::Vector) = v_to_gv(v.(T),T)
v_to_gv(v::LinearInterpolation,T::Vector) = v_to_gv(v.(T),T)

# for Gauss u(t) = 1/(√(σ)*π^(1/4)) * exp( -(t-τ)^2 / (2*σ^2) ) * exp(-i*Δ*t)
function u_to_gu_Gauss(u, τ, σ)
    # u = mode function (Gauss or Gauss*exp(i*ω*t))
    # τ = time delay # σ = width
    ∫u_2_t(t_) = 0.5*(erf( (t_-τ)/σ ) + erf(τ/σ))
    f(t_) = u(t_)' / max(√( 1 - ∫u_2_t(t_) ), _tol_div) # TODO: benchmark
    return f
end
function v_to_gv_Gauss(v, τ, σ)
    # v = mode function (Gauss or Gauss*exp(i*ω*t))
    # τ = time delay # σ = width
    ∫v_2_t(t_) = 0.5*(erf( (t_-τ)/σ ) + erf(τ/σ))
    f(t_) = v(t_)' / max(√( ∫v_2_t(t_) ), _tol_div) # TODO: benchmark
    return f
end

# return g_v(t) from v(t) for multiple output modes
function vi_to_v_i_im1(v_fcts, gv_fcts, T_ls, i)
    @assert i>1
    function multiple_outputs_α(dα,α,p,t) # only for i>1
        for j=1:i-1
            dα[j] = -gv_fcts[j](t)*( v_fcts[i](t) + 
                sum( (gv_fcts[k](t))'*α[k] for k=1:j-1; init=0.0) ) - 
                0.5*abs2(gv_fcts[j](t))*α[j]
        end
    end
    u0 = zeros(ComplexF64, i-1)
    tspan = (T_ls[1], T_ls[end])
    prob = ODEProblem(multiple_outputs_α, u0, tspan)
    sol_α = solve(prob) #; abstol, reltol) # TODO: kwarg
    v_i_im1(t) = v_fcts[i](t) + sum( (gv_fcts[k](t))'*sol_α(t)[k] for k=1:i-1)
    return v_i_im1
end
vi_to_v_i_im1(v_fcts, T_ls, i) = vi_to_v_i_im1(v_fcts, [v_to_gv(v_,T_ls) for v_ in v_fcts], T_ls, i)

# return g_u(t) from u(t) for multiple output modes
function ui_to_u_i_im1(u_fcts, gu_fcts, T_ls, i)
    @assert i>1
    function multiple_inputs_α(dα,α,p,t) # only for i>1
        for j=1:i-1
            dα[j] = -gu_fcts[j](t)*( u_fcts[i](t) - 
                sum( (gu_fcts[i](t))'*α[k] for k=1:j-1; init=0.0) ) + 
                0.5*abs2(gu_fcts[j](t))*α[j]
        end
        # 2 Typos in PRA2020-Kiilerich Eq.(A15): last term "-" → "+" and |g_ui|² → |g_uj|²  
    end
    u0 = zeros(ComplexF64, i-1)
    tspan = (T_ls[1], T_ls[end])
    prob = ODEProblem(multiple_inputs_α, u0, tspan)
    sol_α = solve(prob, Tsit5()) #; abstol, reltol) # TODO: kwarg
    u_i_im1(t) = u_fcts[i](t) - sum( (gu_fcts[k](t))'*sol_α(t)[k] for k=1:i-1)
    return u_i_im1
end
ui_to_u_i_im1(u_fcts, T_ls, i) = ui_to_u_i_im1(u_fcts, [u_to_gu(u_,T_ls) for u_ in u_fcts], T_ls, i)
# TODO: rename, better method (Victor)

# creates matrix of the two time correlation function
function two_time_corr_matrix(T_ls::Vector, ρt::Vector, f::Function, Ls::Function; abstol=1e-6, reltol=1e-6) #kwargs...) TODO
    l_T_ls = length(T_ls)
    @assert l_T_ls == length(ρt)
    Ls_ls = Ls.(T_ls)
    Ls_ls_dag = dagger.(Ls_ls)
    ρ0_ = [Ls_ls[i]*ρt[i] for i=1:l_T_ls] 

    g1_m = zeros(ComplexF64, l_T_ls, l_T_ls)
    for it = 1:l_T_ls-1
        τ_, ρ_bar_τ = timeevolution.master_dynamic(T_ls[it:end], ρ0_[it], f; abstol=abstol, reltol=reltol) #kwargs...) TODO
    
        g1 = [expect(Ls_ls_dag[it+i-1], ρ_bar_τ[i]) for i=1:length(τ_)]
        g1_m[it,it:end] = g1
        g1_m[it:end,it] = adjoint.(g1)
    end
    return g1_m
end
# two_time_corr_matrix for time-independent problems
function two_time_corr_matrix(T_ls::Vector, ρt::Vector, H, J::Vector, Ls; kwargs...)
    l_T_ls = length(T_ls)
    @assert l_T_ls == length(ρt)
    Ls_dag = dagger(Ls)
    ρ0_ = [Ls*ρt[i] for i=1:l_T_ls] 

    g1_m = zeros(ComplexF64, l_T_ls, l_T_ls)
    for it = 1:l_T_ls-1
        τ_, ρ_bar_τ = timeevolution.master(T_ls[it:end], ρ0_[it], H, J; kwargs...) 
    
        g1 = [expect(Ls_dag, ρ_bar_τ[i]) for i=1:length(τ_)]
        g1_m[it,it:end] = g1
        g1_m[it:end,it] = adjoint.(g1)
    end
    return g1_m
end



###########################
### interaction picture ###
###########################

### get matrix M(t)
function get_Mt(A::Function, T)
    T0 = T[1]; Tend = T[end]
    l_A = size(A(T0))[1]
    M0 = diagm(ones(ComplexF64, l_A))
    f_M(u, p, t) = A(t) * u
    prob_M = ODEProblem(f_M, M0, (T0, Tend))
    return solve(prob_M; saveat=T)#, abstol, reltol) # TODO: kwargs
end

# ### get matrix A(t) (for 2x2, 3x3 and 4x4 matrices)
function get_At_2(gu, gv) # assumes coupling u-v
    A(t_) = 0.5*[0 gu(t_)*gv(t_)';   
        -gu(t_)'*gv(t_) 0]
    return A 
end
get_At_2(g_ls) = get_At_2(g_ls...)
function get_At_3(g_ls::Vector) # assumes coupling u-s-v
    @assert length(g_ls) == 3 
    g1 = g_ls[1]; g2 = g_ls[2]; g3 = g_ls[3]
    A(t_) = 0.5*[0 g1(t_)*g2(t_)' g1(t_)*g3(t_)'; 
        -g1(t_)'*g2(t_) 0 g2(t_)*g3(t_)';  
        -g1(t_)'*g3(t_) -g2(t_)'*g3(t_) 0] 
    return A 
end
# gu_ls: list of gu functions for the corresponding operators au
# gv_ls: list of gv functions for the corresponding operators av
function get_At_4(gu_ls::Vector, gv_ls::Vector) 
    @assert length(gu_ls) + length(gv_ls) == 4 
    gu1 = gu_ls[1]; gu2 = gu_ls[2]; gv1 = gv_ls[1]; gv2 = gv_ls[2]
    A(t_) = 0.5*[0 -gu1(t_)*gu2(t_)' gu1(t_)*gv1(t_)' gu1(t_)*gv2(t_)'; 
        gu1(t_)'*gu2(t_) 0 gu2(t_)*gv1(t_)' gu2(t_)*gv2(t_)';  
        -gu1(t_)'*gv1(t_) -gu2(t_)'*gv1(t_) 0 gv1(t_)*gv2(t_)';
        -gu1(t_)'*gv2(t_) -gu2(t_)'*gv2(t_) -gv1(t_)'*gv2(t_) 0]
    return A
end
