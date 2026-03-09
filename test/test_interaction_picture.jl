using QuantumInputOutput
using SecondQuantizedAlgebra
using QuantumOptics
using SymbolicUtils
using LinearAlgebra
using Test

@testset "interaction_picture" begin
    # Parameters
    γ = 1.0
    n_ph = 8

    # Gaussian pulse u(t)
    τ = 1 / γ
    t_p = 4 / γ
    u(t) = 1 / (sqrt(τ) * π^(1/4)) * exp(-0.5 * ((t - t_p) / τ)^2)

    T_end = 12.0
    T = [0:0.002:1;] * T_end

    # Symbolic SLH setup (u -> s -> v)
    hu = FockSpace(:u)
    hs = NLevelSpace(:s, 2)
    hv = FockSpace(:v)
    h = hu ⊗ hs ⊗ hv

    au_sym = Destroy(h, :a_u, 1)
    av_sym = Destroy(h, :a_v, 3)
    σ_sym = Transition(h, :σ, 1, 2, 2)

    gu_sym, γ_sym, gv_sym = rnumbers("gu γ gv")

    G_u = SLH(1, gu_sym * au_sym, 0)
    G_s = SLH(1, sqrt(γ_sym) * σ_sym, 0)
    G_v = SLH(1, gv_sym * av_sym, 0)
    G_cas = ▷(G_u, G_s, G_v)

    H = get_hamiltonian(G_cas)
    L = get_lindblad(G_cas)[1]
    H_uv = get_hamiltonian(▷(G_u, G_v))
    H_int_sym_ = simplify(H - H_uv)

    # Interaction-picture operator substitution
    M(i, j) = cnumber("M_{$(i)$(j)}")
    a0_ls = [au_sym, av_sym]
    la = length(a0_ls)
    a_int_ls = [sum(M(i, j) * a0_ls[j] for j = 1:la) for i = 1:la]
    # int_dict = Dict([a0_ls; adjoint.(a0_ls)] .=> [a_int_ls; adjoint.(a_int_ls)])
    int_dict = Dict(a0_ls .=> a_int_ls)

    H_int_sym = simplify(substitute_operators(H_int_sym_, int_dict))
    L_int_sym = simplify(substitute_operators(L, int_dict))

    # Virtual-cavity couplings
    gu_t = u_to_gu(u, T)
    gv_t = v_to_gv(u, T)

    # Interaction-picture coefficient matrices
    A_uv = interaction_picture_A_2modes(gu_t, gv_t)
    M_num = interaction_picture_M(A_uv, T)
    M_ana = interaction_picture_M_2modes_equal(u, T)

    @test abs(maximum([maximum(abs.(M_num(t))) for t in T]) - 1) < 1e-4
    max_M_err = maximum([maximum(abs.(M_num(t) - M_ana(t))) for t in T])
    @test max_M_err < 5e-4

    # Numerical basis and operators
    bu = FockBasis(n_ph)
    ba = NLevelBasis(2)
    bv = FockBasis(n_ph)
    b = bu ⊗ ba ⊗ bv

    au = destroy(bu) ⊗ one(ba) ⊗ one(bv)
    av = one(bu) ⊗ one(ba) ⊗ destroy(bv)
    σee = one(bu) ⊗ transition(ba, 2, 2) ⊗ one(bv)

    dict_p = Dict(γ_sym => γ)
    M_ls = [M(i, j) for i = 1:la for j = 1:la]
    M_t_ls = [t -> M_num(t)[i, j] for i = 1:la for j = 1:la]
    # M_t_c_ls = [t -> conj(M_num(t)[i, j]) for i = 1:la for j = 1:la]
    # p_t_sym = [gu_sym, gv_sym, M_ls..., conj.(M_ls)...]
    # p_t_num = [gu_t, gv_t, M_t_ls..., M_t_c_ls...]
    p_t_sym = [gu_sym, gv_sym, M_ls...]
    p_t_num = [gu_t, gv_t, M_t_ls...]
    dict_p_t = Dict(p_t_sym .=> p_t_num)

    H_int_QO = translate_qo(H_int_sym, b; parameter = dict_p, time_parameter = dict_p_t)
    L_QO = translate_qo(L_int_sym, b; parameter = dict_p, time_parameter = dict_p_t)

    function input_output_I(t, ρ)
        Ht = H_int_QO(t)
        J = [L_QO(t)]
        return Ht, J, dagger.(J)
    end

    ψ0 = fockstate(bu, n_ph) ⊗ nlevelstate(ba, 1) ⊗ fockstate(bv, 0)
    _, ρt_I = timeevolution.master_dynamic(T, ψ0, input_output_I)

    n_u_t = real.(expect(au' * au, ρt_I))
    n_v_t = real.(expect(av' * av, ρt_I))
    σ_ee = real.(expect(σee, ρt_I))

    σ_ee[end]
    @test abs(σ_ee[end]) < 1e-2
    @test maximum(n_u_t) ≤ n_u_t[1] + 1e-6

    # Non-interaction picture evolution (cascaded u -> s -> v)
    dict_p_s = Dict(γ_sym => γ)
    dict_p_t_s = Dict(gu_sym => gu_t, gv_sym => gv_t)

    H_QO = translate_qo(H, b; parameter = dict_p_s, time_parameter = dict_p_t_s)
    L_QO_S = translate_qo(L, b; parameter = dict_p_s, time_parameter = dict_p_t_s)

    function input_output_S(t, ρ)
        Ht = H_QO(t)
        J = [L_QO_S(t)]
        return Ht, J, dagger.(J)
    end

    _, ρt_S = timeevolution.master_dynamic(T, ψ0, input_output_S)
    σ_ee_S = real.(expect(σee, ρt_S))

    @test maximum(abs.(σ_ee .- σ_ee_S)) < 1e-4

    # Small gamma limit
    γ_small = 1e-4
    dict_p_small = Dict(γ_sym => γ_small)
    H_int_QO_small =
        translate_qo(H_int_sym, b; parameter = dict_p_small, time_parameter = dict_p_t)
    L_QO_small =
        translate_qo(L_int_sym, b; parameter = dict_p_small, time_parameter = dict_p_t)

    function input_output_I_small(t, ρ)
        Ht = H_int_QO_small(t)
        J = [L_QO_small(t)]
        return Ht, J, dagger.(J)
    end

    _, ρt_I_small = timeevolution.master_dynamic(T, ψ0, input_output_I_small)
    n_u_small = real.(expect(au' * au, ρt_I_small))
    n_v_small = real.(expect(av' * av, ρt_I_small))

    @test maximum(n_v_small) < 1e-6
    @test maximum(abs.(n_u_small .- n_u_small[1])) < 5e-3
end
