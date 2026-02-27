using QuantumInputOutput
using SecondQuantizedAlgebra
using QuantumOptics
using LinearAlgebra
using Test

@testset "example_cavity_scattering_minimal" begin
    hu1 = FockSpace(:u1)
    hc1 = FockSpace(:c1)
    hv1 = FockSpace(:v1)
    h = hu1 ⊗ hc1 ⊗ hv1

    au = Destroy(h, :a_u, 1)
    c = Destroy(h, :c, 2)
    av = Destroy(h, :a_v, 3)

    gu, Δ, γ = rnumbers("g_u Δ γ")
    gv = cnumber("g_v")

    G_u = SLH(1, gu * au, 0)
    G_c = SLH(1, √(γ) * c, Δ * c' * c)
    G_v = SLH(1, gv * av, 0)
    G_cas = ▷(G_u, G_c, G_v)

    H = get_hamiltonian(G_cas)
    L = get_lindblad(G_cas)[1]

    γ_ = 1.0
    Δ_ = 0.0
    p_sym = [γ, Δ, gv]
    p_num = [γ_, Δ_, 0.0]
    dict_p = Dict(p_sym .=> p_num)

    T_p = 1/γ_
    T_end = 12T_p
    σ = sqrt(0.5)*T_p
    u1(t) = sqrt(1/(σ*√(2π))*exp( -0.5*(t - 4T_p)^2/σ^2 ))
    T = [0:0.004:1;]*T_end
    ΔT = T[2] - T[1]

    gu_int = u_to_gu(u1, T)
    gu_t(t) = gu_int(t)
    dict_p_t = Dict(gu => gu_t)

    bu1 = FockBasis(2)
    bc1 = FockBasis(2)
    bv1 = FockBasis(2)
    b = bu1 ⊗ bc1 ⊗ bv1

    H_QO = translate(H, b; parameter=dict_p, time_parameter=dict_p_t)
    L_QO = translate(L, b; parameter=dict_p, time_parameter=dict_p_t)

    function input_output_1(t, ρ)
        Ht = H_QO(t)
        J = [L_QO(t)]
        return Ht, J, dagger.(J)
    end

    ψ0 = fockstate(bu1, 1) ⊗ fockstate(bc1, 0) ⊗ fockstate(bv1, 0)
    t_, ρt = timeevolution.master_dynamic(T, ψ0, input_output_1)

    au_qo = translate(au, b)
    c_qo = translate(c, b)
    av_qo = translate(av, b)

    n_u_t = real.(expect(au_qo' * au_qo, ρt))
    n_c_t = real.(expect(c_qo' * c_qo, ρt))

    @test n_u_t[end] < 1e-2
    @test abs(sum(n_u_t .* abs2.(gu_t.(T)))*ΔT - 1) < 1e-2

    # correlation matrix
    Ls(t) = gu_t(t)*au_qo + √(γ_)*c_qo
    g1_m = two_time_corr_matrix(T, ρt, input_output_1, Ls)

    F = eigen(g1_m)
    n_avg = round.(real.(F.values)*ΔT; digits=3)
    modes = F.vectors
    v_mode = (modes[:,end]) / sqrt(ΔT)
    
    @test abs(n_avg[end] - 1) < 1e-2
    @test abs(n_avg[end-1]) < 1e-2

    # output mode
    p_sym_2 = [γ , Δ ]
    p_num_2 = [γ_, Δ_]
    dict_p_2 = Dict(p_sym_2 .=> p_num_2);
    
    # time-depedent coupling for the output mode $v(t)$
    gv_t_ = v_to_gv(v_mode, T)
    gv_t(t) = gv_t_(t)
    gvc_t(t) = conj(gv_t_(t))
    
    dict_p_t_2 = Dict([gu, gv, conj(gv)] .=> [gu_t, gv_t, gvc_t]);
    
    H_QO_2 = translate(H, b; parameter=dict_p_2, time_parameter=dict_p_t_2)
    L_QO_2 = translate(L, b; parameter=dict_p_2, time_parameter=dict_p_t_2)
    function input_output_2(t,ρ)
        H = H_QO_2(t)
        J = [L_QO_2(t)]
        return H, J, dagger.(J)
    end;
    
    # time evolution for the system including the output cavity
    t_2, ρt_2 = timeevolution.master_dynamic(T, ψ0, input_output_2)
    
    n_v1_t = real.(expect(av_qo'*av_qo, ρt_2))
    @test abs(n_v1_t[end] - 1) < 1e-2
    @test abs(sum(n_v1_t .* abs2.(gv_t.(T)))*ΔT - 1) < 1e-2
end   
