using QuantumInputOutput
using SecondQuantizedAlgebra
using QuantumOptics
using QuantumOpticsBase: dagger
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

    @variables Δ::Real γ::Real
    @variables gu::Complex gv::Complex

    G_u = SLH(1, gu' * au, 0)
    G_c = SLH(1, √(γ) * c, Δ * c' * c)
    G_v = SLH(1, gv' * av, 0)
    G_cas = ▷(G_u, G_c, G_v)

    H = hamiltonian(G_cas)
    L = lindblad(G_cas)[1]

    γ_ = 1.0
    Δ_ = 0.0
    p_sym = [γ, Δ, gv]
    p_num = [γ_, Δ_, 0.0]
    dict_p = Dict(p_sym .=> p_num)

    T_p = 1/γ_
    T_end = 12T_p
    σ = sqrt(0.5)*T_p
    u1(t) = 1/(sqrt(σ)*π^(1/4)) * exp(-(t - 4σ)^2 / (2*σ^2))
    T = [0:0.004:1;]*T_end
    ΔT = T[2] - T[1]

    gu_t = coupling_input(u1, T)
    dict_p_t = Dict(gu => gu_t)

    bu1 = FockBasis(2)
    bc1 = FockBasis(2)
    bv1 = FockBasis(2)
    b = bu1 ⊗ bc1 ⊗ bv1

    H_QO = to_numeric(H, b; parameter = dict_p, time_parameter = dict_p_t)
    L_QO = to_numeric(L, b; parameter = dict_p, time_parameter = dict_p_t)

    function input_output_1(t, ρ)
        Ht = H_QO(t)
        J = [L_QO(t)]
        return Ht, J, dagger.(J)
    end

    ψ0 = fockstate(bu1, 1) ⊗ fockstate(bc1, 0) ⊗ fockstate(bv1, 0)
    t_, ρt = timeevolution.master_dynamic(T, ψ0, input_output_1)

    au_qo = to_numeric(au, b)
    c_qo = to_numeric(c, b)
    av_qo = to_numeric(av, b)

    n_u_t = real.(expect(au_qo' * au_qo, ρt))
    n_c_t = real.(expect(c_qo' * c_qo, ρt))

    @test n_u_t[end] < 1e-2
    @test abs(sum(n_u_t .* abs2.(gu_t.(T)))*ΔT - 1) < 1e-2

    # correlation matrix
    Ls(t) = gu_t(t)*au_qo + √(γ_)*c_qo
    g1_m = correlation_matrix(T, ρt, input_output_1, Ls)

    F = eigen(g1_m)
    n_avg = round.(real.(F.values)*ΔT; digits = 3)
    modes = F.vectors
    v_mode = (modes[:, end]) / sqrt(ΔT)

    @test abs(n_avg[end] - 1) < 1e-2
    @test abs(n_avg[end-1]) < 1e-2

    # output mode
    p_sym_2 = [γ, Δ]
    p_num_2 = [γ_, Δ_]
    dict_p_2 = Dict(p_sym_2 .=> p_num_2);

    gv_t = coupling_output(v_mode, T)

    dict_p_t_2 = Dict([gu, gv] .=> [gu_t, gv_t]);

    H_QO_2 = to_numeric(H, b; parameter = dict_p_2, time_parameter = dict_p_t_2)
    L_QO_2 = to_numeric(L, b; parameter = dict_p_2, time_parameter = dict_p_t_2)
    function input_output_2(t, ρ)
        H = H_QO_2(t)
        J = [L_QO_2(t)]
        return H, J, dagger.(J)
    end;

    t_2, ρt_2 = timeevolution.master_dynamic(T, ψ0, input_output_2)

    n_v1_t = real.(expect(av_qo'*av_qo, ρt_2))
    @test abs(n_v1_t[end] - 1) < 1e-2
    @test abs(sum(n_v1_t .* abs2.(gv_t.(T)))*ΔT - 1) < 1e-2
end
