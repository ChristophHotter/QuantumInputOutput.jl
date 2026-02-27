include("/home/christoph/MEGA2/JULIA_FUNCTIONS/QuantumInputOutput/src/QuantumInputOutput.jl")

using PyPlot; pygui(true)
using JLD2
cd(@__DIR__)
abstol = 1e-10
reltol = 1e-10;
pref = "01-1_"

# symbolic Hilbert spaces and operators
N = 3
ha(i) = NLevelSpace("a$(i)", 2)
h = tensor([ha(i) for i=1:N]...)

σ(α,i,j) = Transition(h,"σ_$(α)",i,j,α);

# symbolic paramter
γR(i) = rnumber("γ^{($(i))}_R")
γL(i) = rnumber("γ^{($(i))}_L")
Δ(i) = rnumber("Δ_{$(i)}")
ϕ(i,j) = rnumber("ϕ_{$(i)$(j)}") # phase between QD-i and QD-j
# β(i) = rnumber("β_{$(i)}")
Ein = rnumber("E_{in}")

γR_ = [γR(i) for i=1:N]
γL_ = [γL(i) for i=1:N]
Δ_ = [Δ(i) for i=1:N]
ϕ_ = [ϕ(1,2), ϕ(2,3)]
# β_ = [β(i) for i=1:N]


G_d = SLH(1, Ein, 0)
G_ϕ(i,j) = SLH(exp(1im*ϕ(i,j)), 0, 0)
G_R(i) = SLH(1, √(γR(i))*σ(i,1,2), -Δ(i)*σ(i,2,2))
G_L(i) = SLH(1, √(γL(i))*σ(i,1,2), 0) # don't double-count H0

G_R_t = G_d ▷ G_R(1) ▷ G_ϕ(1,2) ▷ G_R(2) ▷ G_ϕ(2,3) ▷ G_R(3);
G_L_t = G_L(3) ▷ G_ϕ(2,3) ▷ G_L(2) ▷ G_ϕ(1,2) ▷ G_L(1);

G_t = G_R_t ⊞ G_L_t;

S = get_scattering(G_t);
L = get_lindblad(G_t);
H = get_hamiltonian(G_t);


##########################################
### numerical parameters and functions ###
##########################################
# exp: QD1 -> QD3 -> QD2

#####################
nQDs = 2
# transm = true 
######################

for nQDs in 1:3

γn = [0.345, 0.365, 0.388]*2π
# new QD (#2) values are chosen randomly
nQDs==3 && (βn = [0.85, 0.85, 0.95])
nQDs==2 && (βn = [0.85, 0, 0.95]) # to be consistent with the 2QDs simulations
nQDs==1 && (βn = [0.85, 0, 0])

νn = [0.09, 0.09, 0.09]*2π #*0 #TODO
ϕn = [1.2, 1.6]*π 
σSDn = [0.22, 0.26, 0.3]*2π # TODO: monte carlo

Δ_ref = 1 # TODO
Δn = [-0.2, -0.25, -0.3]*2π *Δ_ref 
Ω1n = 0.2 

γRn = γn.*βn/2
γLn = γn.*βn/2
γ_add = γn.*(-βn .+ 1)

En = Ω1n/√(γRn[1])

p_sym = [ γR_; γL_; Δ_; ϕ_; Ein]
p_num = [ γRn; γLn; Δn; ϕn; En]
dict_p = Dict(p_sym .=> p_num)

# numeric basis
ba = NLevelBasis(2)
b = tensor([ba for i=1:N]...)

H_QO = translate(H, b; parameter=dict_p)
L_R_QO = translate(L[1], b; parameter=dict_p)
L_L_QO = translate(L[2], b; parameter=dict_p)
s_QO(α,i,j) =  translate(σ(α,i,j), b; parameter=dict_p)
#
J_add_γ = [√(γ_add[i])*s_QO(i,1,2) for i=1:N]; 
J_add_ν = [√(νn[i]/2)*(s_QO(i,2,2) - s_QO(i,1,1)) for i=1:N]
J_QO = [L_R_QO, L_L_QO, J_add_γ..., J_add_ν...]

Tend = 20/Ω1n
T = [0:0.001:1;]*Tend

for transm in [true, false]
    trre = transm ? "transmission" : "reflection"

    transm ? (L0_QO = L_R_QO) : (L0_QO = L_L_QO)
    L0_QO_dag = dagger(L0_QO)

    ψ0 = tensor([nlevelstate(ba,1) for i=1:N]...)
    t, ρt = timeevolution.master(T, ψ0, H_QO, J_QO; abstol, reltol)
    ρ_ss = ρt[end]


    # close("timeevolution")
    # figure("timeevolution")
    # subplot(211)
    # title("QDs=$(nQDs)")
    # plot(t, real(expect(L0_QO_dag*L0_QO, ρt)))
    # ylabel("⟨L0⁺L0⟩")
    # #
    # subplot(212)
    # plot(t, real(expect(s_QO(1,2,2), ρt)))
    # plot(t, real(expect(s_QO(2,2,2), ρt)), ls="--")
    # ylabel("⟨σ²²⟩")
    # xlabel("t")
    # #
    # tight_layout()
    # savefig("plots/time-evolution/$(pref)_time-evolution.png")

    # correlation function
    T_corr = [0:0.0001:1;]*10
    t_corr, ρt_corr = timeevolution.master(T_corr, L0_QO*ρ_ss*L0_QO_dag, H_QO, J_QO; abstol, reltol) # τ

    G1 = real(expect( L0_QO_dag*L0_QO, ρ_ss))
    Gn_τ(n) = real(expect( L0_QO_dag^(n-1) * L0_QO^(n-1), ρt_corr))
    gn_τ(n) = Gn_τ(n)/G1^n

    t_corr_plot = [-reverse(t_corr); t_corr]
    gn_τ_plot(n) = [reverse(gn_τ(n)); gn_τ(n)]
    Gn_τ_plot(n) = [reverse(Gn_τ(n)); Gn_τ(n)]

    n_max = 4
    cd(@__DIR__)
    close("g(n)τ")
    figure("g(n)τ", figsize=(4,8))
    for n=2:n_max
        subplot(n_max-1,1,n-1)
        n==2 && title("$(trre); QDs=$(nQDs)")
        plot(t_corr_plot, gn_τ_plot(n))
        grid(true)
        ylabel("g($(n))")
        xlim(-5,5)

        # save data
        gn_τ_plot_n = gn_τ_plot(n)
        Gn_τ_plot_n = Gn_τ_plot(n)
        @save "data/$(pref)_QDs=$(nQDs)_$(trre)_Delta-ref=$(Δ_ref)__gn_τ_plot_n.jld2" gn_τ_plot_n
        @save "data/$(pref)_QDs=$(nQDs)_$(trre)_Delta-ref=$(Δ_ref)__Gn_τ_plot_n.jld2" Gn_τ_plot_n
    end
    xlabel("τ")
    tight_layout()
    savefig("plots/$(pref)_QDs=$(nQDs)_$(trre)_Delta-ref=$(Δ_ref)_gn-tau.png")
end # transm

end # nQDs