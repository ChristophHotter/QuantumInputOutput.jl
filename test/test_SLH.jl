using QuantumInputOutput
using SecondQuantizedAlgebra
using QuantumOptics
using QuantumOpticsBase: dagger
using SymbolicUtils
using FunctionWrappers: FunctionWrapper
using StaticArrays
using Test

@testset "SLH" begin
    hu1 = FockSpace(:u1)
    hc1 = FockSpace(:c1)
    hv1 = FockSpace(:v1)
    h = hu1 ⊗ hc1 ⊗ hv1

    au = Destroy(h, :a_u, 1)
    c = Destroy(h, :c, 2)
    av = Destroy(h, :a_v, 3)

    @variables gu::Real Δ::Real γ::Real
    @variables gv::Complex

    G_u = SLH(1, gu'*au, 0) # input cavity
    G_c = SLH(1, √(γ)*c, Δ*c'c) # system cavity
    G_v = SLH(1, gv'*av, 0) # output cavity

    G_c_S = scattering(G_c)
    G_c_L = jump_operator(G_c)
    G_c_H = hamiltonian(G_c)

    @test G_c_S == G_c.scattering
    @test G_c_L == G_c.jump_operator
    @test G_c_H == G_c.hamiltonian
    @test_deprecated lindblad(G_c)

    @test G_c_S isa SMatrix{1,1}
    @test G_c_L isa SVector{1}

    SLH(1, [√(γ)*c], Δ*c'c)
    @test isequal(G_c, SLH(1, [√(γ)*c], Δ*c'c))

    @testset "simple_cascade" begin
        G1 = G_u ▷ G_c
        @test scattering(G1) isa SMatrix{1,1}
        @test iszero(simplify(jump_operator(G1)[1] - simplify(gu'*au + √(γ)*c)))
        expected_H = simplify(
            hamiltonian(G_c) - 1im/2*((√(γ)*c)'*(1)*gu'*au - (gu'*au)'*(1)*(√(γ)*c)),
        )
        @test iszero(simplify(hamiltonian(G1) - expected_H))

        G2 = cascade(G_u, G_c, G_v)
        G3 = G_u ▷ G_c ▷ G_v

        @test isequal(G2, G1 ▷ G_v)
        @test isequal(G2, ▷(G1, G_v))
        @test isequal(G2, G3)

        @test iszero(simplify(jump_operator(G2)[1] - (gu'*au + √(γ)*c + gv'*av)))
    end

    @testset "simple_concatenate" begin
        G1 = SLH(1, gu'*au, 0)
        G2 = SLH(1, √(γ)*c, Δ*c'c)

        Gc = concatenate(G1, G2)

        @test scattering(Gc) isa SMatrix{2,2}
        @test size(scattering(Gc)) == (2, 2)
        @test length(jump_operator(Gc)) == 2
        @test isequal(jump_operator(Gc)[1], gu'*au)
        @test isequal(jump_operator(Gc)[2], √(γ)*c)
        @test isequal(hamiltonian(Gc), Δ*c'c)

        Gc2 = G1 ⊞ G2
        @test isequal(Gc, Gc2)
    end

    @testset "two-port symbolic cascade" begin
        hu2 = FockSpace(:u2)
        hv2 = FockSpace(:v2)
        h2 = hu1 ⊗ hu2 ⊗ hv1 ⊗ hv2

        au1 = Destroy(h2, :a_u1, 1)
        au2 = Destroy(h2, :a_u2, 2)
        av1 = Destroy(h2, :a_v1, 3)
        av2 = Destroy(h2, :a_v2, 4)

        @variables gu1::Real gu2::Real gv1::Real gv2::Real t::Real r::Real

        G_bs = SLH([r t; t -r], [0, 0], 0)
        G_out = SLH(1, gv1' * av1, 0) ⊞ SLH(1, gv2' * av2, 0)
        G_cas = G_bs ▷ G_out

        @test length(jump_operator(G_cas)) == 2
        @test iszero(simplify(jump_operator(G_cas)[1] - (gv1' * av1)))
        @test iszero(simplify(jump_operator(G_cas)[2] - (gv2' * av2)))
    end

    @testset "numeric type stability" begin
        bc = FockBasis(4)
        a_op = destroy(bc)
        H_s = sparse(0.5 * dagger(a_op) * a_op)
        L_s = sparse(sqrt(1.0) * a_op)
        gu_f(t) = exp(-t^2) * sparse(a_op)
        gv_f(t) = exp(-(t - 2)^2) * sparse(a_op)

        @testset "static SLH concreteness" begin
            G = SLH(1, L_s, H_s)
            @test eltype(jump_operator(G)) === typeof(L_s)
            @test typeof(hamiltonian(G)) === typeof(H_s)
        end

        @testset "time-dep SLH wraps into FunctionWrapper" begin
            G_td = SLH(1, gu_f, H_s)
            @test eltype(jump_operator(G_td)) <: FunctionWrapper
            @test typeof(hamiltonian(G_td)) <: FunctionWrapper
        end

        @testset "cascade preserves FunctionWrapper" begin
            G1 = SLH(1, gu_f, H_s)
            G2 = SLH(1, gv_f, H_s)
            G_cas = G1 ▷ G2
            @test eltype(jump_operator(G_cas)) <: FunctionWrapper
            @test typeof(hamiltonian(G_cas)) <: FunctionWrapper
            @inferred jump_operator(G_cas)[1](0.5)
        end

        @testset "cascade mixed static/time-dep wraps uniformly" begin
            G_cas = SLH(1, L_s, H_s) ▷ SLH(1, gu_f, H_s)
            @test eltype(jump_operator(G_cas)) <: FunctionWrapper
            @test eltype(jump_operator(G_cas)) !== Any
        end

        @testset "concatenation mixed static/time-dep wraps uniformly" begin
            G_cat = SLH(1, L_s, H_s) ⊞ SLH(1, gu_f, H_s)
            LT = eltype(jump_operator(G_cat))
            @test LT <: FunctionWrapper
            @test LT !== Any
            @inferred jump_operator(G_cat)[1](0.5)
            @inferred jump_operator(G_cat)[2](0.5)
        end

        @testset "concatenation static stays static" begin
            G_cat = SLH(1, L_s, H_s) ⊞ SLH(1, L_s, H_s)
            @test eltype(jump_operator(G_cat)) === typeof(L_s)
            @test !(eltype(jump_operator(G_cat)) <: FunctionWrapper)
        end

        @testset "FunctionWrapper call is inferred" begin
            G_td = SLH(1, gu_f, H_s)
            l = jump_operator(G_td)[1]
            @inferred l(0.5)
        end

        @testset "_op_type extracts type from FunctionWrapper SLH" begin
            G_td = SLH(1, gu_f, H_s)
            @test QuantumInputOutput._op_type(G_td) === typeof(H_s)
        end

        @testset "_op_type returns nothing for static SLH" begin
            G_s = SLH(1, L_s, H_s)
            @test QuantumInputOutput._op_type(G_s) === nothing
        end

        @testset "SLH with only plain closures errors" begin
            bare_f(t) = t * ones(ComplexF64, 5, 5)
            bare_g(t) = (1 - t) * ones(ComplexF64, 5, 5)
            @test_throws ErrorException SLH([1 0; 0 1], [bare_f, bare_g], bare_g)
        end

        @testset "feedback preserves FunctionWrapper type" begin
            G1 = SLH(1, gu_f, H_s)
            G2 = SLH(1, gv_f, H_s)
            G_cat = G1 ⊞ G2
            G_fb = feedback(G_cat, 2, 1)
            @test eltype(jump_operator(G_fb)) <: FunctionWrapper
            @test typeof(hamiltonian(G_fb)) <: FunctionWrapper
        end
    end
end
