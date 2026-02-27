using QuantumInputOutput
using QuantumCumulants
using SymbolicUtils
using QuantumOpticsBase
using Test

@testset "translate" begin
    @rnumbers κ_L κ_R Δ g γ
    @cnumbers E
    Natoms = 2

    hc = FockSpace(:cavity)
    ha_ = NLevelSpace("a",2)
    h = hc ⊗ ha_

    a = Destroy(h,:a,1) # cavity 
    σ(i,j) = Transition(h,"σ",i,j,2) # two-level atom 

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
    dict_p_t2 = Dict( [E, conj(E)] .=> [E_t, E_t_c] )

    @testset "kwarg_operators" begin
        @test isequal(translate(Δ, bc1; parameter=dict_p1, operators=ops_dict), one(bc1)*Δn)
        @test isequal(translate(2, bc1; parameter=dict_p1, operators=ops_dict), one(bc1)*2)
        @test isequal(translate(a, bc1; parameter=dict_p1, operators=ops_dict), a_QO)
        @test isequal(translate(a*3, bc1; parameter=dict_p1, operators=ops_dict), a_QO*3)
        F1 = translate(a*3, bc1; parameter=dict_p1, time_parameter=dict_p_t2, operators=ops_dict)
        @test isa(F1, Function)
        @test isequal(F1(0.1),a_QO*3)
        F2 = translate(a*E, bc1; parameter=dict_p1, time_parameter=dict_p_t2, operators=ops_dict)
        @test isequal(F2(0.1),a_QO*E_t(0.1))
        F2(0.1)
        a_QO*E_t(0.1)
        @test isequal(translate(Δ*a'a, bc1; parameter=dict_p1, operators=ops_dict), Δn*dagger(a_QO)*a_QO)
    end

    bc1 = FockBasis(4)
    ba = NLevelBasis(2)
    b = bc1 ⊗ ba

    a_QO2 = destroy(bc1) ⊗ one(ba)
    σ_QO(i,j) = to_numeric(σ(i,j),b)
    
    @test isequal(a_QO2, dense(to_numeric(a,b)))
    @test isequal(translate(Δ, b; parameter=dict_p1), one(b)*Δn)
    F3 = translate(Δ, b; parameter=dict_p1, time_parameter=dict_p_t2)
    @test sum(abs.((F3(4) - one(b)*Δn).data)) < 1e-8
    F4 = translate(a*3*conj(E) + Δ*σ(2,2), b; parameter=dict_p1, time_parameter=dict_p_t2)
    @test sum(abs.((F4(0.2) - dense(a_QO2*3*E_t_c(0.2) + Δn*σ_QO(2,2))).data)) < 1e-8
    F5 = translate(a*conj(E) + Δ*σ(2,2), b; parameter=dict_p1, time_parameter=dict_p_t2)
    @test sum(abs.((F5(0.2) - dense(a_QO2*E_t_c(0.2) + Δn*σ_QO(2,2))).data)) < 1e-8
    F5_ = translate(a*conj(E), b; parameter=dict_p1, time_parameter=dict_p_t2)
    @test sum(abs.((F5_(0.2) - dense(a_QO2*E_t_c(0.2))).data)) < 1e-8
    F6 = translate(conj(E), b; parameter=dict_p1, time_parameter=dict_p_t2)
    @test sum(abs.((F6(0.2) - dense(E_t_c(0.2)*one(b))).data)) < 1e-8
    F7 = translate(conj(E) + Δ*σ(2,2), b; parameter=dict_p1, time_parameter=dict_p_t2)
    @test sum(abs.((F7(0.2) - dense(E_t_c(0.2)*one(b) + Δn*σ_QO(2,2))).data)) < 1e-8
    @test_throws MethodError translate(conj(E), b; parameter=dict_p1)
    # F8 = translate(E^2, b; parameter=dict_p1, time_parameter=dict_p_t2) # TODO

    @testset "time_parameter_normalization" begin
        dict_p_t_num = Dict([E] .=> [2.5])
        F_num = translate(a*E, b; parameter=dict_p1, time_parameter=dict_p_t_num)
        @test sum(abs.((F_num(0.2) - dense(a_QO2 * 2.5)).data)) < 1e-8
    end
end

