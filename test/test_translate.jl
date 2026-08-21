using QuantumInputOutput
using SecondQuantizedAlgebra
using SymbolicUtils
using QuantumOpticsBase
using QuantumOpticsBase: dagger, static_operator, TimeDependentSum
using Test

mat(F, t) = dense(static_operator(F(t)))

@testset "translate" begin
    @variables κ_L::Real κ_R::Real Δ::Real g::Real γ::Real
    @variables E::Complex E1::Complex E2::Complex
    Natoms = 2

    hc = FockSpace(:cavity)
    ha_ = NLevelSpace(:a, 2)
    h = hc ⊗ ha_

    a = Destroy(h, :a, 1) # cavity
    σ(i, j) = Transition(h, :σ, i, j, 2) # two-level atom

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
    dict_p_t2 = Dict(E => E_t)


    @testset "kwarg_operators" begin
        @test isequal(
            to_numeric(Δ, bc1; parameter = dict_p1, operators = ops_dict),
            one(bc1)*Δn,
        )
        @test isequal(
            to_numeric(2, bc1; parameter = dict_p1, operators = ops_dict),
            one(bc1)*2,
        )
        @test isequal(to_numeric(a, bc1; parameter = dict_p1, operators = ops_dict), a_QO)
        @test isequal(
            to_numeric(a*3, bc1; parameter = dict_p1, operators = ops_dict),
            a_QO*3,
        )
        F1 = to_numeric(
            a*3,
            bc1;
            parameter = dict_p1,
            time_parameter = dict_p_t2,
            operators = ops_dict,
        )
        @test F1 isa TimeDependentSum
        @test sum(abs.((mat(F1, 0.1) - dense(a_QO*3)).data)) < 1e-12
        F2 = to_numeric(
            a*E,
            bc1;
            parameter = dict_p1,
            time_parameter = dict_p_t2,
            operators = ops_dict,
        )
        @test sum(abs.((mat(F2, 0.1) - dense(a_QO*E_t(0.1))).data)) < 1e-12
        # association of the scalar prefactor differs from `Δn*A'*A`, compare numerically
        @test sum(
            abs.(
                (
                    dense(
                        to_numeric(Δ*a'a, bc1; parameter = dict_p1, operators = ops_dict),
                    ) - dense(Δn*dagger(a_QO)*a_QO)
                ).data,
            ),
        ) < 1e-12
    end

    bc1 = FockBasis(4)
    ba = NLevelBasis(2)
    b = bc1 ⊗ ba

    a_QO2 = destroy(bc1) ⊗ one(ba)
    σ_QO(i, j) = to_numeric(σ(i, j), b)

    @test isequal(
        a_QO2'σ_QO(1, 2),
        dense(
            to_numeric(
                a'σ(1, 2),
                b;
                operators = Dict([a, σ(1, 2)] .=> [a_QO2, σ_QO(1, 2)]),
                adjoint_ops = true,
            ),
        ),
    )
    @test isequal(a_QO2, dense(to_numeric(a, b)))
    @test sum(
        abs.((dense(to_numeric(Δ, b; parameter = dict_p1)) - dense(one(b)*Δn)).data),
    ) < 1e-12
    F3 = to_numeric(Δ, b; parameter = dict_p1, time_parameter = dict_p_t2)
    @test sum(abs.((mat(F3, 4) - dense(one(b)*Δn)).data)) < 1e-8
    F4 = to_numeric(
        a*3*conj(E) + Δ*σ(2, 2),
        b;
        parameter = dict_p1,
        time_parameter = dict_p_t2,
    )
    @test sum(abs.((mat(F4, 0.2) - dense(a_QO2*3*E_t_c(0.2) + Δn*σ_QO(2, 2))).data)) < 1e-8
    F5 = to_numeric(
        a*conj(E) + Δ*σ(2, 2),
        b;
        parameter = dict_p1,
        time_parameter = dict_p_t2,
    )
    @test sum(abs.((mat(F5, 0.2) - dense(a_QO2*E_t_c(0.2) + Δn*σ_QO(2, 2))).data)) < 1e-8
    F5_ = to_numeric(a*conj(E), b; parameter = dict_p1, time_parameter = dict_p_t2)
    @test sum(abs.((mat(F5_, 0.2) - dense(a_QO2*E_t_c(0.2))).data)) < 1e-8
    F6 = to_numeric(conj(E), b; parameter = dict_p1, time_parameter = dict_p_t2)
    @test sum(abs.((mat(F6, 0.2) - dense(E_t_c(0.2)*one(b))).data)) < 1e-8
    F7 = to_numeric(conj(E) + Δ*σ(2, 2), b; parameter = dict_p1, time_parameter = dict_p_t2)
    @test sum(abs.((mat(F7, 0.2) - dense(E_t_c(0.2)*one(b) + Δn*σ_QO(2, 2))).data)) < 1e-8
    # a bare symbolic scalar without a numeric/time value cannot be translated
    @test_throws ArgumentError to_numeric(conj(E), b; parameter = dict_p1)
    F8 = to_numeric(E^2, b; parameter = dict_p1, time_parameter = dict_p_t2)
    @test sum(abs.((mat(F8, 0.2) - dense(E_t(0.2)^2*one(b))).data)) < 1e-8


    @testset "time_parameter_normalization" begin
        dict_p_t_num = Dict([E] .=> [2.5])
        F_num = to_numeric(a*E, b; parameter = dict_p1, time_parameter = dict_p_t_num)
        @test sum(abs.((mat(F_num, 0.2) - dense(a_QO2 * 2.5)).data)) < 1e-8
    end

    @testset "multiple_time_prefactors" begin
        E1_t(t) = 1.2 + 0.3im + 0.5t
        E2_t(t) = 0.7 - 0.1im + 0.2t
        E2_t_c(t) = conj(E2_t(t))
        dict_p_t_multi = Dict([E1, conj(E2)] .=> [E1_t, E2_t_c])
        F_multi = to_numeric(
            a * E1 * conj(E2),
            b;
            parameter = dict_p1,
            time_parameter = dict_p_t_multi,
        )
        expected = dense(a_QO2 * E1_t(0.4) * E2_t_c(0.4))
        @test sum(abs.((mat(F_multi, 0.4) - expected).data)) < 1e-8
    end

    @testset "to_numeric vector overload" begin
        ops = [a, Δ, a * E]
        translated_ops = to_numeric(
            ops,
            bc1;
            parameter = dict_p1,
            time_parameter = dict_p_t2,
            operators = ops_dict,
        )
        @test length(translated_ops) == length(ops)
        @test all(op -> op isa TimeDependentSum, translated_ops)
        @test sum(abs.((mat(translated_ops[1], 0.3) - dense(a_QO)).data)) < 1e-12
        @test sum(abs.((mat(translated_ops[2], 0.3) - dense(one(bc1) * Δn)).data)) < 1e-12
        @test sum(abs.((mat(translated_ops[3], 0.3) - dense(a_QO * E_t(0.3))).data)) < 1e-12
    end

    @testset "to_numeric SLH overload" begin
        G_sym = SLH(1, [sqrt(κ_R) * a * E], Δ * a' * a)
        H_QO, J_QO = to_numeric(
            G_sym,
            bc1;
            parameter = dict_p1,
            time_parameter = dict_p_t2,
            operators = ops_dict,
        )

        @test H_QO isa TimeDependentSum
        @test length(J_QO) == 1
        @test J_QO[1] isa TimeDependentSum
        @test sum(abs.((mat(H_QO, 0.4) - dense(Δn * dagger(a_QO) * a_QO)).data)) < 1e-8
        @test sum(abs.((mat(J_QO[1], 0.4) - dense(sqrt(κ_Rn) * a_QO * E_t(0.4))).data)) <
              1e-8
    end

    @testset "substitute operators qmul" begin
        h2 = FockSpace(:h2)
        a = Destroy(h2, :a)
        @variables c1::Complex c2::Complex c3::Complex
        a_1 = 2 * c2 * a + c3*a
        a_2 = c2 * a * a
        dict_sub = Dict(a => a_1)
        dict_sub2 = Dict(a => a_2)

        x = a*a*c1 + a*c3
        y = substitute(x, dict_sub)
        y2 = substitute(x, dict_sub2)

        @test iszero(y - (a_1*a_1*c1 + a_1*c3))
        @test iszero(y2 - (a_2*a_2*c1 + a_2*c3))

        @test substitute(5, dict_sub) == 5
    end

    @testset "substitute adjoint in args_nc" begin
        h2 = FockSpace(:h2)
        a = Destroy(h2, :a)
        b = Destroy(h2, :b)
        ad = a'
        bd = b'
        @variables g::Complex

        # product with adjoint operator: g * a† * a
        op = g * ad * a
        dict_sub = Dict(a => b)

        result = substitute(op, dict_sub)
        expected = g * bd * b
        @test iszero(result - expected)
    end

    @testset "complex-valued time function for a real-typed parameter" begin
        # A coefficient `im*gR` (with `gR` declared `::Real`) fed a complex-valued time
        # function used to throw `MethodError: complex(::Int64, ::ComplexF64)` because
        # `build_function` emitted a `complex(re, im)` call. The compile path now combines
        # the real and imaginary parts with `im`, so it works for complex-valued inputs.
        h3 = FockSpace(:h3)
        a3 = Destroy(h3, :a)
        b3 = FockBasis(3)
        a3_QO = destroy(b3)
        @variables gR::Real
        gR_t(t) = 1.0 + 2.0im
        F = to_numeric(im * gR * a3, b3; time_parameter = Dict(gR => gR_t))
        @test F isa TimeDependentSum
        @test sum(abs.((mat(F, 0.0) - dense((im * (1.0 + 2.0im)) * a3_QO)).data)) < 1e-8
    end

    @testset "PulseCoupling as time_parameter" begin
        # A `PulseCoupling` from a coupling constructor must plug into `to_numeric`'s
        # `time_parameter` without the conflicting-arity error a raw interpolation triggers.
        hcp = FockSpace(:cp)
        acp = Destroy(hcp, :a)
        bcp = FockBasis(3)
        acp_QO = destroy(bcp)
        @variables gu::Complex
        Tg = collect(0.0:0.05:1.0)
        umode(t) = exp(-(t - 0.5)^2)
        gu_t = coupling_input(umode, Tg)
        @test gu_t isa PulseCoupling
        Fpc = to_numeric(
            acp * gu,
            bcp;
            time_parameter = Dict(gu => gu_t),
            operators = Dict([acp, acp'] .=> [acp_QO, dagger(acp_QO)]),
        )
        @test Fpc isa TimeDependentSum
        @test sum(abs.((mat(Fpc, 0.3) - dense(gu_t(0.3) * acp_QO)).data)) < 1e-10
    end
end
