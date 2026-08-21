using QuantumInputOutput
using SecondQuantizedAlgebra
using QuantumOptics
using QuantumOpticsBase: dagger
using SymbolicUtils
using Symbolics: Symbolics
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

    @variables gu_sym::Complex gv_sym::Complex
    @variables γ_sym::Real

    G_u = SLH(1, gu_sym' * au_sym, 0)
    G_s = SLH(1, sqrt(γ_sym) * σ_sym, 0)
    G_v = SLH(1, gv_sym' * av_sym, 0)
    G_cas = ▷(G_u, G_s, G_v)

    H = hamiltonian(G_cas)
    Jop = jump_operator(G_cas)[1]
    H_uv = hamiltonian(▷(G_u, G_v))
    H_int_sym_ = simplify(H - H_uv)

    # Interaction-picture operator substitution
    M(i, j) = Symbolics.variable(Symbol("M_{$(i)$(j)}"); T = Complex{Real})
    a0_ls = [au_sym, av_sym]
    la = length(a0_ls)
    a_int_ls = [sum(M(i, j) * a0_ls[j] for j = 1:la) for i = 1:la]
    int_dict = Dict(a0_ls .=> a_int_ls)

    H_int_sym = simplify(substitute(H_int_sym_, int_dict))
    J_int_sym = simplify(substitute(Jop, int_dict))

    # Virtual-cavity couplings
    gu_t = coupling_input(u, T)
    gv_t = coupling_output(u, T)

    # Interaction-picture coefficient matrices
    A_uv = coupling_matrix((gu_t, gv_t))
    sol_M = solve_mode_evolution(A_uv, T)
    M_num = t -> sol_M(t)
    M_ana = solve_mode_evolution_symmetric(u, T)

    @test abs(maximum([maximum(abs.(M_num(t))) for t in T]) - 1) < 1e-4
    max_M_err = maximum([maximum(abs.(M_num(t) - M_ana(t))) for t in T])
    @test max_M_err < 5e-4

    # Complex couplings: A_ij = g_i g_j^* / 2 for i < j, -g_j^* g_i / 2 for i > j
    let g1 = 1 + 2im, g2 = 3 + 4im, g3 = -1 + 0.5im
        A2 = coupling_matrix((g1, g2))(0.0)
        @test A2 ≈ 0.5 * [0 g1*conj(g2); -conj(g1)*g2 0]

        A3 = coupling_matrix(g1, g2, g3)(0.0)
        A3_expected =
            0.5 * [
                0 g1*conj(g2) g1*conj(g3)
                -conj(g1)*g2 0 g2*conj(g3)
                -conj(g1)*g3 -conj(g2)*g3 0
            ]
        @test A3 ≈ A3_expected
        @test A3' ≈ -A3  # anti-Hermitian generator
    end

    @static if VERSION > v"1.12.0"
        # coupling_matrix type stability
        @test A_uv(0.5) isa SMatrix{2,2,ComplexF64}
        @inferred A_uv(0.5)
        A_uv(0.0)  # warmup
        @test (@allocated A_uv(0.5)) == 0

        # also for a tuple of mixed element type (constant next to interpolant)
        A_mixed = coupling_matrix(gu_t, 1.0 + 0.5im, gv_t)
        @test A_mixed(0.5) isa SMatrix{3,3,ComplexF64}
        @inferred A_mixed(0.5)
        A_mixed(0.0)  # warmup
        @test (@allocated A_mixed(0.5)) == 0
    end

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
    p_t_sym = [gu_sym, gv_sym, M_ls...]
    p_t_num = [gu_t, gv_t, M_t_ls...]
    dict_p_t = Dict(p_t_sym .=> p_t_num)

    H_int_QO = to_numeric(H_int_sym, b; parameter = dict_p, time_parameter = dict_p_t)
    J_QO = to_numeric(J_int_sym, b; parameter = dict_p, time_parameter = dict_p_t)

    function input_output_I(t, ρ)
        Ht = H_int_QO(t)
        J = [J_QO(t)]
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

    H_QO = to_numeric(H, b; parameter = dict_p_s, time_parameter = dict_p_t_s)
    J_QO_S = to_numeric(Jop, b; parameter = dict_p_s, time_parameter = dict_p_t_s)

    function input_output_S(t, ρ)
        Ht = H_QO(t)
        J = [J_QO_S(t)]
        return Ht, J, dagger.(J)
    end

    _, ρt_S = timeevolution.master_dynamic(T, ψ0, input_output_S)
    σ_ee_S = real.(expect(σee, ρt_S))

    @test maximum(abs.(σ_ee .- σ_ee_S)) < 1e-4

    # Small gamma limit
    γ_small = 1e-4
    dict_p_small = Dict(γ_sym => γ_small)
    H_int_QO_small =
        to_numeric(H_int_sym, b; parameter = dict_p_small, time_parameter = dict_p_t)
    J_QO_small =
        to_numeric(J_int_sym, b; parameter = dict_p_small, time_parameter = dict_p_t)

    function input_output_I_small(t, ρ)
        Ht = H_int_QO_small(t)
        J = [J_QO_small(t)]
        return Ht, J, dagger.(J)
    end

    _, ρt_I_small = timeevolution.master_dynamic(T, ψ0, input_output_I_small)
    n_u_small = real.(expect(au' * au, ρt_I_small))
    n_v_small = real.(expect(av' * av, ρt_I_small))

    @test maximum(n_v_small) < 1e-6
    @test maximum(abs.(n_u_small .- n_u_small[1])) < 5e-3
end
