using QuantumInputOutput
using SecondQuantizedAlgebra
using SymbolicUtils
using StaticArrays
using Test

@testset "SLH" begin
    hu1 = FockSpace(:u1)
    hc1 = FockSpace(:c1)
    hv1 = FockSpace(:v1)
    h = hu1 ⊗ hc1 ⊗ hv1

    au = Destroy(h, :a_u, 1)
    c = Destroy(h, :c, 2)
    av = Destroy(h, :a_v, 3);

    gu, Δ, γ = rnumbers("g_u Δ, γ")
    gv = cnumber("g_v");

    G_u = SLH(1, gu'*au, 0) # input cavity
    G_c = SLH(1, √(γ)*c, Δ*c'c) # system cavity
    G_v = SLH(1, gv'*av, 0) # output cavity

    G_c_S = scattering(G_c)
    G_c_L = lindblad(G_c)
    G_c_H = hamiltonian(G_c)

    @test G_c_S == G_c.scattering
    @test G_c_L == G_c.lindblad
    @test G_c_H == G_c.hamiltonian

    @test G_c_S isa SMatrix{1,1}
    @test G_c_L isa SVector{1}

    SLH(1, [√(γ)*c], Δ*c'c)
    @test isequal(G_c, SLH(1, [√(γ)*c], Δ*c'c))

    @testset "simple_cascade" begin
        G1 = G_u ▷ G_c
        @test scattering(G1) isa SMatrix{1,1}
        @test isequal(lindblad(G1)[1], simplify(gu'*au + √(γ)*c))
        expected_H = simplify(hamiltonian(G_c) - 1im/2*((√(γ)*c)'*(1)*gu'*au - (gu'*au)'*(1)*(√(γ)*c)))
        @test iszero(simplify(hamiltonian(G1) - expected_H))

        G2 = cascade(G_u, G_c, G_v)
        G3 = G_u ▷ G_c ▷ G_v

        @test isequal(G2, G1 ▷ G_v)
        @test isequal(G2, ▷(G1, G_v))
        @test isequal(G2, G3)

        @test iszero(simplify(lindblad(G2)[1] - (gu'*au + √(γ)*c + gv'*av)))
    end

    @testset "simple_concatenate" begin
        G1 = SLH(1, gu'*au, 0)
        G2 = SLH(1, √(γ)*c, Δ*c'c)

        Gc = concatenate(G1, G2)

        @test scattering(Gc) isa SMatrix{2,2}
        @test size(scattering(Gc)) == (2, 2)
        @test length(lindblad(Gc)) == 2
        @test isequal(lindblad(Gc)[1], gu'*au)
        @test isequal(lindblad(Gc)[2], √(γ)*c)
        @test isequal(hamiltonian(Gc), Δ*c'c)

        Gc2 = G1 ⊞ G2
        @test isequal(Gc, Gc2)
    end

end
