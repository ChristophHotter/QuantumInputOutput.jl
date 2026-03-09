using QuantumInputOutput
using SecondQuantizedAlgebra
using SymbolicUtils
using QuantumOpticsBase
using Test

@testset "translate" begin
    @rnumbers κ_L κ_R Δ g γ
    @cnumbers E E1 E2
    Natoms = 2

    hc = FockSpace(:cavity)
    ha_ = NLevelSpace("a", 2)
    h = hc ⊗ ha_

    a = Destroy(h, :a, 1) # cavity 
    σ(i, j) = Transition(h, "σ", i, j, 2) # two-level atom 

    bc1 = FockBasis(4)
    a_QO = destroy(bc1)
    ops_dict = Dict([a, a'] .=> [a_QO, dagger(a_QO)])

    En = 0.5
    κ_Rn = 1.5
    κ_Ln = 1.0
    Δn = 0.2

    p_sym = [κ_R, κ_L, Δ]
    p_num = [κ_Rn, κ_Ln, Δn]
    dict_p1 = Dict(p_sym .=> p_num)

    E_t(t) = 2*t + 1im
    E_t_c(t) = conj(E_t(t))
    # dict_p_t2 = Dict( [E, conj(E)] .=> [E_t, E_t_c] )
    dict_p_t2 = Dict(E => E_t)


    @testset "kwarg_operators" begin
        @test isequal(
            translate(Δ, bc1; parameter = dict_p1, operators = ops_dict),
            one(bc1)*Δn,
        )
        @test isequal(
            translate(2, bc1; parameter = dict_p1, operators = ops_dict),
            one(bc1)*2,
        )
        @test isequal(translate(a, bc1; parameter = dict_p1, operators = ops_dict), a_QO)
        @test isequal(
            translate(a*3, bc1; parameter = dict_p1, operators = ops_dict),
            a_QO*3,
        )
        F1 = translate(
            a*3,
            bc1;
            parameter = dict_p1,
            time_parameter = dict_p_t2,
            operators = ops_dict,
        )
        @test isa(F1, Function)
        @test isequal(F1(0.1), a_QO*3)
        F2 = translate(
            a*E,
            bc1;
            parameter = dict_p1,
            time_parameter = dict_p_t2,
            operators = ops_dict,
        )
        @test isequal(F2(0.1), a_QO*E_t(0.1))
        F2(0.1)
        a_QO*E_t(0.1)
        @test isequal(
            translate(Δ*a'a, bc1; parameter = dict_p1, operators = ops_dict),
            Δn*dagger(a_QO)*a_QO,
        )
    end

    bc1 = FockBasis(4)
    ba = NLevelBasis(2)
    b = bc1 ⊗ ba

    a_QO2 = destroy(bc1) ⊗ one(ba)
    σ_QO(i, j) = to_numeric(σ(i, j), b)

    @test isequal(a_QO2, dense(to_numeric(a, b)))
    @test isequal(translate(Δ, b; parameter = dict_p1), one(b)*Δn)
    F3 = translate(Δ, b; parameter = dict_p1, time_parameter = dict_p_t2)
    @test sum(abs.((F3(4) - one(b)*Δn).data)) < 1e-8
    F4 = translate(
        a*3*conj(E) + Δ*σ(2, 2),
        b;
        parameter = dict_p1,
        time_parameter = dict_p_t2,
    )
    @test sum(abs.((F4(0.2) - dense(a_QO2*3*E_t_c(0.2) + Δn*σ_QO(2, 2))).data)) < 1e-8
    F5 =
        translate(a*conj(E) + Δ*σ(2, 2), b; parameter = dict_p1, time_parameter = dict_p_t2)
    @test sum(abs.((F5(0.2) - dense(a_QO2*E_t_c(0.2) + Δn*σ_QO(2, 2))).data)) < 1e-8
    F5_ = translate(a*conj(E), b; parameter = dict_p1, time_parameter = dict_p_t2)
    @test sum(abs.((F5_(0.2) - dense(a_QO2*E_t_c(0.2))).data)) < 1e-8
    F6 = translate(conj(E), b; parameter = dict_p1, time_parameter = dict_p_t2)
    @test sum(abs.((F6(0.2) - dense(E_t_c(0.2)*one(b))).data)) < 1e-8
    F7 = translate(conj(E) + Δ*σ(2, 2), b; parameter = dict_p1, time_parameter = dict_p_t2)
    @test sum(abs.((F7(0.2) - dense(E_t_c(0.2)*one(b) + Δn*σ_QO(2, 2))).data)) < 1e-8
    @test_throws MethodError translate(conj(E), b; parameter = dict_p1)
    F8 = translate(E^2, b; parameter = dict_p1, time_parameter = dict_p_t2)
    @test sum(abs.((F8(0.2) - dense(E_t(0.2)^2*one(b))).data)) < 1e-8


    @testset "time_parameter_normalization" begin
        dict_p_t_num = Dict([E] .=> [2.5])
        F_num = translate(a*E, b; parameter = dict_p1, time_parameter = dict_p_t_num)
        @test sum(abs.((F_num(0.2) - dense(a_QO2 * 2.5)).data)) < 1e-8
    end

    @testset "multiple_time_prefactors" begin
        E1_t(t) = 1.2 + 0.3im + 0.5t
        E2_t(t) = 0.7 - 0.1im + 0.2t
        E2_t_c(t) = conj(E2_t(t))
        dict_p_t_multi = Dict([E1, conj(E2)] .=> [E1_t, E2_t_c])
        F_multi = translate(
            a * E1 * conj(E2),
            b;
            parameter = dict_p1,
            time_parameter = dict_p_t_multi,
        )
        expected = dense(a_QO2 * E1_t(0.4) * E2_t_c(0.4))
        @test sum(abs.((F_multi(0.4) - expected).data)) < 1e-8
    end

    @testset "substitute_operators_qmul" begin
        h2 = FockSpace(:h2)
        a = Destroy(h2, :a, 1)
        @cnumbers c1 c2 c3
        a_1 = 2 * c2 * a + c3*a
        a_2 = c2 * a * a
        dict_sub = Dict(a => a_1)
        dict_sub2 = Dict(a => a_2)

        x = a*a*c1 + a*c3
        y = substitute_operators(x, dict_sub)
        y2 = substitute_operators(x, dict_sub2)

        @test isequal(simplify(y - (a_1*a_1*c1 + a_1*c3)), 0)
        @test isequal(simplify(y2 - (a_2*a_2*c1 + a_2*c3)), 0)
    end
end
