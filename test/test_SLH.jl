using QuantumInputOutput
using QuantumCumulants
using SymbolicUtils
using Test

# @testset "SLH" begin
    hu1 = FockSpace(:u1)
    hc1 = FockSpace(:c1)
    hv1 = FockSpace(:v1)
    h = hu1 ⊗ hc1 ⊗ hv1

    au = Destroy(h,:a_u,1)
    c = Destroy(h,:c,2)
    av = Destroy(h,:a_v,3);

    gu, Δ, γ = rnumbers("g_u Δ, γ")
    gv = cnumber("g_v");

    G_u = SLH(1, gu*au, 0) # input cavity 
    G_c = SLH(1, √(γ)*c, Δ*c'c) # system cavity
    G_v = SLH(1, gv*av, 0) # output cavity

    G_c_S = get_scattering(G_c)
    G_c_L = get_lindblad(G_c)
    G_c_H = get_hamiltonian(G_c)
    
    @test G_c_S == G_c.scattering
    @test G_c_L == G_c.lindblad
    @test G_c_H == G_c.hamiltonian

    SLH(1, [√(γ)*c], Δ*c'c)
    @test isequal(G_c, SLH(1, [√(γ)*c], Δ*c'c))
    @test isequal(G_c, SLH(ones(1,1), [√(γ)*c], Δ*c'c))

    @test isequal(ones(1,1), G_c_S)
    @test isequal([√(γ)*c], G_c_L)
    @test isequal(Δ*c'c, G_c_H)

    @testset "simple_cascade" begin
        G1 = G_u ▷ G_c
        @test isequal(G1.scattering, ones(1,1))
        @test isequal(G1.lindblad[1], simplify(gu*au + √(γ)*c))
        @test isequal(G1.hamiltonian, simplify(G_c_H - 1im/2*( (√(γ)*c)'*(1)*gu*au - (gu*au)'*(1)*(√(γ)*c) )))
        
        G2 = cascade(G_u, G_c, G_v)
        
        @test isequal(G2, G1 ▷ G_v)
        @test isequal(G2, ▷(G1,G_v))
        @test iszero(simplify(G2.lindblad[1] - (gu*au + √(γ)*c + gv*av)))
    end

    # @testset "simple_concatenation" begin
    #     #    
    # end


# end

    1

    # @test first(typeof(G).parameters) == Complex{Real}
    # @test first(typeof(ω).parameters) == Real

    # @test SymbolicUtils.getmetadata(ω, Symbolics.VariableSource) == (:RealParameter, :ω)
    # @test SymbolicUtils.getmetadata(G, Symbolics.VariableSource) == (:Parameter, :G)

    # @testset "rnumbers_cnumbers" begin
    #     @test isequal(rnumbers(:r1, :r2), (r1, r2))
    #     @test isequal(rnumbers("r1"), (r1,))
    #     @test isequal(rnumber("r1"), (r1))
    #     @test isequal(rnumber(:r1), (r1))
    #     @test isequal(cnumber("c"), cnumber(:c))
    # end

    # @testset "one_zero" begin
    #     @test one(ω) == 1.0
    #     @test one(G) == 1.0 + 0.0im
    #     @test zero(ω) == 0.0
    #     @test zero(G) == 0.0 + 0.0im
    # end

    # @testset "adjoint" begin
    #     @test isequal(adjoint(ω), ω)
    #     @test isequal(adjoint(3ω), 3ω)
    #     @test isequal(adjoint(G), conj(G))
    #     @test_broken isequal(conj(3*G), 3*conj(G))
    #     @test isequal(SecondQuantizedAlgebra._conj(3*G), 3*conj(G))
    #     @test isequal(simplify(adjoint(G*G1)), conj(G*G1))
    #     @test isequal(SecondQuantizedAlgebra._adjoint(exp(1im*ϕ)*r1), exp(-1im*ϕ)*r1)
    #     @test isequal(SecondQuantizedAlgebra._conj(exp(1im*ϕ)*r1), exp(-1im*ϕ)*r1)
    #     @test simplify(SecondQuantizedAlgebra._adjoint(exp(1im*ϕ)*exp(-1im*ϕ))) == 1.0
    # end
# end
