### #TODO:

# # TODO: same in other file: delete conj fcts
# # TODO: delete nv term in other file! not correct in Δn!


# # Input-Output Analysis of Quantum Dot SUPER Excitation
#
# This script collects and demonstrates several approaches used in [PAPER REF].
#
# We provide a quantum pulse description of the SUPER excitation scheme [T. K. Bracht et al., PRX Quantum 2, 040354 (2021)] (https://doi.org/10.1103/PRXQuantum.2.040354)
# and show that SUPER is a three-photon process [Q. W. Richter et al., Physical Review Research 7, 013079 (2025)] (https://doi.org/10.1103/PhysRevResearch.7.013079), 
# [L. Vannucci and N. Gregersen, Optics Express Vol. 32, Issue 20, pp. 35381-35394 (2024)] (https://doi.org/10.1364/OE.533998). 
# Due to large photon numbers, we use QuantumCumulants.jl [D. Plankensteiner et al., Quantum 6, 617 (2022)] (https://doi.org/10.22331/q-2022-01-04-617).
#
# We also propose full master-equation treatments for large-α coherent states in [PAPER REF]; here we show the “Displacement Frame” approach.
#
# Additionally we demonstrate how suitable output modes can be determined utilizing a correlation-based approach
#
# 1) Two-Input-Two-Output (2I-2O) SUPER Cascade
#
# Two Gaussian coherent-state pulses are emitted by two virtual input cavities, drive a TLS, and are collected by matched output cavities (v_i(t) = u_i(t)).
# The difference between input and output photon numbers per pulse reveals the three-photon nature.
#
# 2) Two-Input-Two-Output (2I-2O) SUPER Concatenation
#
# The same physics can be captured by concatenating channels: for each pulse i, input cavity i -> TLS -> output cavity i.
# The SLH description has two inputs and two outputs; photon and system dynamics recover the multi-photon character.
#
# 3) Interaction Picture and Displacement Frame of 2I-2O Cascade
#
# We transform the 2I-2O cascade into the cavity–cavity interaction picture [V. R. Christiansen et al., Physical Review A 107, 013706 (2023)] (https://doi.org/10.1103/PhysRevA.107.013706).
# We then enter a displacement frame that splits coherent backgrounds from quantum fluctuations, reducing the Fock spaces and enabling full master-equation treatment.
#
# 4) Input-Output SUPER Cascade
#
# We find a single coherent mode representing the SUPER superposition. This compressed description hides the explicit multi-photon nature, but greatly speeds up the QuantumCumulants.jl-based evaluation of the two-time correlation matrix g^(1)(t_1, t_2).
# Diagonalizing g^(1) yields the dominant output modes of the Input–TLS SLH triplet and of the full I–O SUPER cascade.
# # TODO: re-order 

using QuantumInputOutput
using SecondQuantizedAlgebra
using QuantumCumulants
using ModelingToolkitBase
using OrdinaryDiffEq
using LinearAlgebra
using QuantumOptics
using NumericalIntegration
using JLD2
using Plots
using LaTeXStrings

abstol = 1e-10
reltol = 1e-10
cd(@__DIR__)
solution_values(sol, op, eqs) = get_solution(sol, op, eqs).(sol.t)
# parameter_map(sys, ps, vals) =
# Dict(p => v for (p, v) in zip(ps, vals) if any(isequal(p), parameters(sys)))

###################################################
############# 1) 2I-2O SUPER Cascade ##############
###################################################

# Hilbert spaces
hu2 = FockSpace(:u2) # virtual input cavity (u2)
hu1 = FockSpace(:u1) # virtual input cavity (u1)
hs1 = NLevelSpace(:atom, 2) # TLS
hv1 = FockSpace(:v1) # virtual output cavity (v1)
hv2 = FockSpace(:v2) # virtual output cavity (v2)
h = tensor(hu2, hu1, hs1, hv1, hv2)

# Operators
au2 = Destroy(h, :a_u2, 1)
au1 = Destroy(h, :a_u1, 2)
s(i, j) = Transition(h, :s, i, j)
av1 = Destroy(h, :a_v1, 4)
av2 = Destroy(h, :a_v2, 5)

# Symbolic parameters: decay rate and virtual cavity couplings
@variables γ::Real gu1::Number gu2::Number gv1::Number gv2::Number
@independent_variables t # Symbolic time variable

###################################################
############# SUPER Pulse Parameters ##############
###################################################

# Time grid
dt = 1e-4;
Tend = 12*2
T = [dt:dt:Tend;];
ΔT = T[2] - T[1]

# Numerical pulse parameters
γ_ = 1e-2 # TLS decay rate
Δ1_ = -2π*1.934 # detuning pulse 1
Δ2_ = -2π*4.634 # detuning pulse 2
α1_ = 22.65π # pulse area 1
α2_ = 19.29π # pulse area 2
σ1_ = 2.4 # temporal FWHM pulse 1
σ2_ = 3.04 # temporal FWHM pulse 2
τ1_ = 0.0 + 12 # time shift pulse 1
τ2_ = -0.73 + 12 # time shift pulse 2

# Classical Gaussian pulses
Ω1(t) = 0.5*1/(√(2π))*α1_/σ1_*exp(-(t-τ1_)^2/(2*σ1_^2))*exp(-1im*Δ1_*t)
Ω2(t) = 0.5*1/(√(2π))*α2_/σ2_*exp(-(t-τ2_)^2/(2*σ2_^2))*exp(-1im*Δ2_*t)

###################################################
##### Classical SUPER Scheme (reference) ##########
###################################################

@register_symbolic Ω(t)
@register_symbolic Ω_c(t)

Ω(t) = Ω1(t) + Ω2(t)
Ω_c(t) = conj(Ω(t))

H_cl = -(Ω_c(t)*s(1, 2) + Ω(t)*s(2, 1))
J_cl = [s(1, 2)];
rates_cl = [γ]

ops_cl = [s(2, 2), s(1, 2)]
eqs_cl = meanfield(ops_cl, H_cl, J_cl; rates = rates_cl, iv = t)
sys_cl = mtkcompile(System(eqs_cl; name = :sys_cl))

u0_cl = zeros(ComplexF64, length(eqs_cl))
u0_p_map_cl = merge(Dict(unknowns(sys_cl) .=> u0_cl), Dict(γ => γ_))
prob_cl = ODEProblem(sys_cl, u0_p_map_cl, (dt, Tend))
sol_cl = solve(prob_cl, Tsit5(); abstol, reltol)
t_cl = sol_cl.t
s22_cl = real.(solution_values(sol_cl, s(2, 2), eqs_cl));

###################################################
####### IOT and QuantumCumulants.jl section #######
###################################################

# SLH triplets
G_u2 = SLH(1, gu2*au2, 0) # input cavity 2
G_u1 = SLH(1, gu1*au1, 0) # input cavity 1
G_2lvl = SLH(1, √(γ)*s(1, 2), 0) # 2-level system
G_v1 = SLH(1, gv1*av1, 0) # output cavity 1
G_v2 = SLH(1, gv2*av2, 0) # output cavity 2

# 2I-2O cascade SLH triplet
G_cas = ▷(G_u2, G_u1, G_2lvl, G_v1, G_v2) # cascade

# Hamiltonian and Lindbladian
Hcas = hamiltonian(G_cas)
Lcas = lindblad(G_cas)[1]
Lcasd = adjoint(Lcas)

# Time-dependent couplings
@register_symbolic gu1_t(t)
@register_symbolic gu2_t(t)
@register_symbolic gv1_t(t)
@register_symbolic gv2_t(t)
# Complex conjugates
@register_symbolic gu1_c_t(t)
@register_symbolic gu2_c_t(t)
@register_symbolic gv1_c_t(t)
@register_symbolic gv2_c_t(t)

g_ls = [gu2, gu1, gv1, gv2]
gt_ls = [gu2_t(t), gu1_t(t), gv1_t(t), gv2_t(t)]
gct_ls = [gu2_c_t(t), gu1_c_t(t), gv1_c_t(t), gv2_c_t(t)]
dict_gt = Dict([g_ls; conj.(g_ls)] .=> [gt_ls; gct_ls]);

# Insert time-dependence
Hcas_t = substitute(Hcas, dict_gt);
Lcas_t = substitute(Lcas, dict_gt);
Lcasd_t = substitute(Lcasd, dict_gt);

# Normalized input modes
u1(t_) = 1/(√(σ1_)*π^(1/4)) * exp(-(t_ - τ1_)^2 / (2*σ1_^2)) * exp(-1im*Δ1_*t_)
u2(t_) = 1/(√(σ2_)*π^(1/4)) * exp(-(t_ - τ2_)^2 / (2*σ2_^2)) * exp(-1im*Δ2_*t_)

# Coupling functions
gu1_t_ = coupling_input(u1, T)
gu1_t(t) = gu1_t_(t)
gv1_t_ = coupling_output(u1, T)
gv1_t(t) = gv1_t_(t)
gu1_c_t(t) = conj(gu1_t(t))
gv1_c_t(t) = conj(gv1_t(t))

# Cascade-modified couplings (effective modes)
u_fcts = [u1, u2]
u2_eff = effective_input_mode(u_fcts, T, 2; abstol, reltol)
gu2_t_ = coupling_input(u2_eff, T)
gu2_t(t) = gu2_t_(t)
gu2_c_t(t) = conj(gu2_t(t));

v_fcts = [u1, u2] # output modes = input modes
v2_eff = effective_output_mode(v_fcts, T, 2; abstol, reltol)
gv2_t_ = coupling_output(v2_eff, T)
gv2_t(t) = gv2_t_(t)
gv2_c_t(t) = conj(gv2_t(t));

# First-order cumulant expansion
order = 1
ops = [au1, au2, s(2, 2), s(2, 1), av1, av2]
eqs = meanfield(ops, Hcas_t, [Lcas_t]; Jdagger = [Lcasd_t], order = order, iv = t)

# Classical-to-coherent amplitude relation
α1_io = α1_ / (2*√(2)*π^(1/4)*√(σ1_*γ_)) # field 1
α2_io = α2_ / (2*√(2)*π^(1/4)*√(σ2_*γ_)) # field 2
u0 = [α1_io, α2_io, 0, 0, 0, 0.0im]

# Solve ODE system
sys = mtkcompile(System(eqs; name = :sys))
u0_p_map_cas = Dict([unknowns(sys); γ] .=> [u0; γ_])
prob_cas = ODEProblem(sys, u0_p_map_cas, (dt, Tend))
sol = solve(prob_cas, Tsit5(); abstol, reltol);

# Expectation values
t_cas = sol.t # time vector
s22_cas = solution_values(sol, s(2, 2), eqs);
nu1_cas = abs2.(solution_values(sol, au1, eqs))
nu2_cas = abs2.(solution_values(sol, au2, eqs))
nv1_cas = abs2.(solution_values(sol, av1, eqs))
nv2_cas = abs2.(solution_values(sol, av2, eqs))

# Plots
plot_font = 18
common = (;
    xlims = (t_cas[1]-0.01, t_cas[end]),
    tickfontsize = plot_font,
    guidefontsize = plot_font,
    legendfontsize = plot_font,
)
p1 = plot(
    t_cl,
    s22_cl;
    color = :blue,
    label = L"\mathrm{classical}",
    ylabel = L"P",
    xticks = ([0, 10, 20], ["", "", ""]),
    yticks = ([0, 0.5, 1], [L"0.0", L"0.5", L"1.0"]),
    legend = :topleft,
    common...,
)
plot!(p1, t_cas, real.(s22_cas); color = :red, ls = :dash, label = L"\mathrm{cascade}")
p2 = plot(
    t_cas,
    nu1_cas;
    color = :blue,
    label = L"\langle \hat{n}_{u_1} \rangle",
    ylabel = L"\langle\hat{n}_i\rangle",
    xticks = ([0, 10, 20], ["", "", ""]),
    yticks = ([0, 5e3, 10e3, 15e3], [L"0", L"5\cdot10^3", L"10\cdot10^3", L"15\cdot10^3"]),
    ylims = (nv1_cas[1]-0.5e3, nu1_cas[1]+0.5e3),
    legend = :left,
    common...,
)
plot!(p2, t_cas, nv1_cas; color = :red, label = L"\langle \hat{n}_{v_1} \rangle")
plot!(
    p2,
    t_cas,
    nu2_cas;
    color = :blue,
    ls = :dash,
    label = L"\langle \hat{n}_{u_2} \rangle",
)
plot!(
    p2,
    t_cas,
    nv2_cas;
    color = :red,
    ls = :dash,
    label = L"\langle \hat{n}_{v_2} \rangle",
)
p3 = plot(
    t_cas,
    nu1_cas .+ nv1_cas .- nu1_cas[1];
    color = :blue,
    label = L"\mathrm{mode~1}",
    ylabel = L"\langle\Delta\hat{n}_i\rangle",
    xlabel = L"\gamma t",
    xticks = ([0, 10, 20], latexstring.([0, 10, 20])),
    yticks = ([-4, -2, 0, 2, 4], latexstring.([-4, -2, 0, 2, 4])),
    legend = :topleft,
    common...,
)
plot!(p3, t_cas, nu2_cas .+ nv2_cas .- nu2_cas[1]; color = :red, label = L"\mathrm{mode~2}")
pl1 = plot(p1, p2, p3; layout = (3, 1), size = (800, 800))
display(pl1)
# savefig(pl1, "plot_2I2O_cascade.pdf")


####################################################################
##### 3) 2I-2O SUPER Cascade interaction picture and Displacement ##
####################################################################

# Interaction picture: cavity dynamics
H_uv = hamiltonian(▷(G_u2, G_u1, G_v1, G_v2))
H_int_sym_ = simplify(Hcas - H_uv)

M(i, j) = Symbolics.variable(Symbol("M_{$(i)$(j)}"); T = Number)
Mc(i, j) = Symbolics.variable(Symbol("Mc_{$(i)$(j)}"); T = Number)

a0_ls = [au2, au1, av1, av2]
la = length(a0_ls)
a_int_ls = [sum(M(i, j)*a0_ls[j] for j = 1:la) for i = 1:la]
# a_c_int_ls = [sum(Mc(i,j)*a0_ls[j]' for j = 1:la) for i = 1:la] # # TODO
a_c_int_ls = [sum(conj(M(i, j))*a0_ls[j]' for j = 1:la) for i = 1:la] # # TODO

int_dict = Dict([a0_ls; adjoint.(a0_ls)] .=> [a_int_ls; a_c_int_ls])

H_int_sym = substitute_operators(H_int_sym_, int_dict)
L_int_sym = simplify(substitute_operators(Lcas, int_dict))
Ld_int_sym = simplify(substitute_operators(Lcasd, int_dict))

# M-matrix
Mat = Matrix{Any}(undef, la, la)
Matc = Matrix{Any}(undef, la, la)
for i = 1:la, j = 1:la
    name = Symbol("Ma_$(i)$(j)")
    namec = Symbol("Mac_$(i)$(j)")

    @eval @register_symbolic $name(t)
    @eval @register_symbolic $namec(t)

    Mat[i, j] = getfield(Main, name)(t)
    Matc[i, j] = getfield(Main, namec)(t)
end
Mat_ls = [Mat[i, j] for i = 1:la for j = 1:la]
Mat_conls = [Matc[i, j] for i = 1:la for j = 1:la]

M_ls = [M(i, j) for i = 1:la for j = 1:la]
Mc_ls = [Mc(i, j) for i = 1:la for j = 1:la]

# Time-evolution matrix M(t)
A_uv = coupling_matrix((gu2_t_, gu1_t_, gv1_t_, gv2_t_))
M_t = solve_mode_evolution(A_uv, T)
Mc_t(t) = adjoint(M_t(t))

for i = 1:la, j = 1:la
    fname = Symbol("Ma_$(i)$(j)")
    fnamec = Symbol("Mac_$(i)$(j)")

    @eval begin
        $fname(t) = M_t(t)[$i, $j]
        $fnamec(t) = Mc_t(t)[$i, $j]
    end
end

dict_Mt = Dict([M_ls..., Mc_ls...] .=> [Mat_ls..., Mat_conls...])
dict_gt_Mt = merge(dict_gt, dict_Mt)

H_int_sym_t = substitute(H_int_sym, dict_gt_Mt)
L_int_sym_t = substitute(L_int_sym, dict_gt_Mt)
Ld_int_sym_t = substitute(Ld_int_sym, dict_gt_Mt)

eqs_int = meanfield(
    ops,
    H_int_sym_t,
    [L_int_sym_t];
    Jdagger = [Ld_int_sym_t],
    order = order,
    iv = t,
);

# Solve ODE system in interaction picture
sys_int = mtkcompile(System(eqs_int; name = :sysI))
u0_p_map_int = Dict([unknowns(sys_int); γ] .=> [u0; γ_])
prob_int = ODEProblem(sys_int, u0_p_map_int, (dt, Tend))
sol_int = solve(prob_int, Tsit5(); abstol, reltol);

# Expectation values
t_int = sol_int.t
s22_int = solution_values(sol_int, s(2, 2), eqs_int);
nu1_int = abs2.(solution_values(sol_int, au1, eqs_int))
nu2_int = abs2.(solution_values(sol_int, au2, eqs_int))
nv1_int = abs2.(solution_values(sol_int, av1, eqs_int))
nv2_int = abs2.(solution_values(sol_int, av2, eqs_int))

### Displacement frame ###

# QuantumOptics.jl Hilbert space for fluctuations # TODO: Why is a cutoff of 1 enough?
bu1 = FockBasis(1) # vacuum fluctuation mode 1
bu2 = FockBasis(1) # vacuum fluctuation mode 2
bs1 = NLevelBasis(2) # TLS
bv1 = FockBasis(1) # vacuum mode 1
bv2 = FockBasis(1) # vacuum mode 2
b = tensor([bu2, bu1, bs1, bv1, bv2]...)

ψu1 = fockstate(bu1, 0)
ψu2 = fockstate(bu2, 0)
ψs1 = nlevelstate(bs1, 1)
ψv1 = fockstate(bv1, 0)
ψv2 = fockstate(bv2, 0)
ψ0 = tensor(ψu2, ψu1, ψs1, ψv1, ψv2)

# Operators
anu2 = destroy(bu2)
anu1 = destroy(bu1)
anv1 = destroy(bv1)
anv2 = destroy(bv2)
s12 = embed(b, 3, transition(bs1, 1, 2))
s22 = embed(b, 3, transition(bs1, 2, 2))

# Mode operator vector
avn = [embed(b, 1, anu2), embed(b, 2, anu1), embed(b, 4, anv1), embed(b, 5, anv2)]

# Coherent amplitude shifts (input modes only)
αvec = [α2_io, α1_io, 0.0, 0.0]

# Displacement-frame Hamiltonian and jumps
function Ham_displaced(t, ρ)
    Mall = M_t(t)
    gall = [gu2_t(t), gu1_t(t), gv1_t(t), gv2_t(t)]

    # Modified coupling vectors for H and L terms
    gn = [gall[1], gall[2], -gall[3], -gall[4]]
    gnc = conj.(gn)
    gl = [gall[1], gall[2], gall[3], gall[4]]
    glc = conj.(gl)

    # Subspaces and local annihilation operators
    B = [bu2, bu1, bv1, bv2]
    Av = [anu2, anu1, anv1, anv2]

    # Displaced operators: a -> b + α (b = fluctuation mode)
    avn_disp = []
    for k = 1:4
        if k < 3
            push!(avn_disp, embed(b, k, (Av[k] + αvec[k]*one(B[k]))))
        else
            push!(avn_disp, embed(b, k+1, (Av[k] + αvec[k]*one(B[k]))))
        end
    end

    # Interaction pieces (IP + displacement)
    A_gn = transpose(gn) * Mall * avn_disp
    A_gnc = avn_disp' * Mall' * gnc

    A_gl = transpose(gl) * Mall * avn_disp
    A_glc = avn_disp' * Mall' * glc

    # Displacement Hamiltonian
    Ht = sqrt(γ_) * 1im/2 * (A_gn * s12' - s12 * A_gnc)

    # Jumps (displacement Lindbladian)
    J = [sqrt(γ_) * s12 + A_gl]
    Jd = [sqrt(γ_) * s12' + A_glc]

    return Ht, J, Jd
end

# Master equation in displacement frame
t_df, ρt_df = timeevolution.master_dynamic(T, ψ0, Ham_displaced; abstol, reltol);

# Expectation values
s22_df = real.(expect(s22, ρt_df))
nu2_df = real.(
    expect(
        embed(
            b,
            1,
            (anu2 + α2_io*identityoperator(bu2))'*(anu2 + α2_io*identityoperator(bu2)),
        ),
        ρt_df,
    ),
)
nu1_df = real.(
    expect(
        embed(
            b,
            2,
            (anu1 + α1_io*identityoperator(bu1))'*(anu1 + α1_io*identityoperator(bu1)),
        ),
        ρt_df,
    ),
)
nv1_df = real.(expect(embed(b, 4, anv1'*anv1), ρt_df))
nv2_df = real.(expect(embed(b, 5, anv2'*anv2), ρt_df))


# Plot net loss/gain in interaction picture and displacement frames
pl4 = plot(
    t_int,
    nu1_int .+ nv1_int .- nu1_int[1];
    color = :blue,
    label = L"\mathrm{mode~1~(int.)}",
    ylabel = L"\langle\Delta\hat{n}_i\rangle",
    xlabel = L"\gamma t",
    xlims = (t_df[1]-0.01, t_df[end]),
    yticks = ([-2, -1, 0, 1], latexstring.([-2, -1, 0, 1])),
    legend = :right,
    tickfontsize = plot_font,
    guidefontsize = plot_font,
    legendfontsize = plot_font,
    size = (800, 400),
)
plot!(
    pl4,
    t_int,
    nu2_int .+ nv2_int .- nu2_int[1];
    color = :red,
    label = L"\mathrm{mode~2~(int.)}",
)
plot!(
    pl4,
    t_df,
    nu1_df .+ nv1_df .- nu1_df[1];
    color = :green,
    ls = :dash,
    label = L"\mathrm{mode~1~(dis.)}",
)
plot!(
    pl4,
    t_df,
    nu2_df .+ nv2_df .- nu2_df[1];
    color = :orange,
    ls = :dash,
    label = L"\mathrm{mode~2~(dis.)}",
)
display(pl4)
# savefig(pl4, "plot_int_dis.pdf")


###################################################
######### 2) 2I-2O SUPER Concatenation ############
###################################################

# TODO: No substitution with γ/2! - interpretation: decay rate twice

# SLH triplets
G_u2 = SLH(1, gu2*au2, 0) # input cavity 2
G_u1 = SLH(1, gu1*au1, 0) # input cavity 1
G_2lvl = SLH(1, √(γ)*s(1, 2), 0) # 2-level system # TODO
G_v1 = SLH(1, gv1*av1, 0) # output cavity 1
G_v2 = SLH(1, gv2*av2, 0) # output cavity 2

# Channel-wise cascades
G_ch1 = ▷(G_u2, G_2lvl, G_v2) # pulse 1
G_ch2 = ▷(G_u1, G_2lvl, G_v1) # pulse 2
Gcon = ⊞(G_ch1, G_ch2) # concatenation

# Hamiltonian and Lindblad operators (two channels)
Hcon = hamiltonian(Gcon)
Lch1g = lindblad(Gcon)[1]
Lch2g = lindblad(Gcon)[2]
Ldch1g = adjoint(Lch1g)
Ldch2g = adjoint(Lch2g)

# Split decay γ -> γ/2 (symmetric into two channels)
# dict_gam = Dict(γ => γ/2); # TODO
dict_gam = Dict(γ => γ);

# Apply decay adjustment
Lch1 = substitute(Lch1g, dict_gam) # TODO
Lch2 = substitute(Lch2g, dict_gam)
Ldch1 = substitute(Ldch1g, dict_gam)
Ldch2 = substitute(Ldch2g, dict_gam)

# Insert time-dependence
Hcon_t = substitute(Hcon, dict_gt);
Lch1_t = substitute(Lch1, dict_gt);
Lch2_t = substitute(Lch2, dict_gt);
Ldch1_t = substitute(Ldch1, dict_gt);
Ldch2_t = substitute(Ldch2, dict_gt);

gu2_t_ = coupling_input(u2, T)
gu2_t(t) = gu2_t_(t)
gv2_t_ = coupling_output(u2, T)
gv2_t(t) = gv2_t_(t)

gu2_c_t(t) = conj(gu2_t(t))
gv2_c_t(t) = conj(gv2_t(t))

# First-order cumulant expansion
order = 1
ops = [au1, au2, s(2, 2), s(2, 1), av1, av2]
eqs_con = meanfield(
    ops,
    Hcon_t,
    [Lch1_t, Lch2_t];
    Jdagger = [Ldch1_t, Ldch2_t],
    order = order,
    iv = t,
)

# Solve ODE system
sys_con = mtkcompile(System(eqs_con; name = :sys_con))

u0_p_map_con = Dict([unknowns(sys_con); γ] .=> [u0; γ_]) # TODO
prob_con = ODEProblem(sys_con, u0_p_map_con, (dt, Tend))
sol_con = solve(prob_con, Tsit5(); abstol, reltol);

# Expectation values
t_con = sol_con.t
s22_con = solution_values(sol_con, s(2, 2), eqs_con)
nu1_con = abs2.(solution_values(sol_con, au1, eqs_con))
nu2_con = abs2.(solution_values(sol_con, au2, eqs_con))
nv1_con = abs2.(solution_values(sol_con, av1, eqs_con))
nv2_con = abs2.(solution_values(sol_con, av2, eqs_con))


common = (;
    xlims = (t_con[1]-0.01, t_con[end]),
    tickfontsize = plot_font,
    guidefontsize = plot_font,
    legendfontsize = plot_font,
)
p1 = plot(
    t_cl,
    real.(s22_cl);
    color = :blue,
    label = L"\mathrm{cascade}",
    ylabel = L"P",
    xticks = ([0, 10, 20], ["", "", ""]),
    legend = :topleft,
    common...,
)
plot!(p1, t_con, real.(s22_con); color = :red, ls = :dash, label = L"\mathrm{concat.}")
p2 = plot(
    t_con,
    nu1_con;
    color = :blue,
    label = L"\langle \hat{n}_{u_1} \rangle",
    ylabel = L"\langle\hat{n}_i\rangle",
    xticks = ([0, 10, 20], ["", "", ""]),
    ylims = (nv1_con[1]-0.5e3, nu1_con[1]+0.5e3),
    legend = :left,
    common...,
)
plot!(p2, t_con, nv1_con; color = :red, label = L"\langle \hat{n}_{v_1} \rangle")
plot!(
    p2,
    t_con,
    nu2_con;
    color = :blue,
    ls = :dash,
    label = L"\langle \hat{n}_{u_2} \rangle",
)
plot!(
    p2,
    t_con,
    nv2_con;
    color = :red,
    ls = :dash,
    label = L"\langle \hat{n}_{v_2} \rangle",
)
p3 = plot(
    t_con,
    nu1_con .+ nv1_con .- nu1_con[1];
    color = :blue,
    label = L"\mathrm{mode~1}",
    ylabel = L"\langle\Delta\hat{n}_i\rangle",
    xlabel = L"\gamma t",
    yticks = ([-2, -1, 0, 1, 2], latexstring.([-2, -1, 0, 1, 2])),
    legend = :bottomleft,
    common...,
)
plot!(p3, t_con, nu2_con .+ nv2_con .- nu2_con[1]; color = :red, label = L"\mathrm{mode~2}")
pl2 = plot(p1, p2, p3; layout = (3, 1), size = (800, 800))
display(pl2)
# savefig(pl2, "plot_concat.pdf")


###################################################
########## 4) I-O SUPER Cascade ###################
###################################################

# Single-mode SUPER description
dt = 1e-5
Tend = 12*2
T = [dt:dt:Tend;]
ΔT = T[2] - T[1]

# Hilbert space
hu = FockSpace(:u) # virtual input cavity 
hv = FockSpace(:v) # virtual output cavity 
hs = tensor(hu, hs1, hv)

# Operators
aun = Destroy(hs, :a_u, 1)
σ(i, j) = Transition(hs, :σ, i, j)
avn = Destroy(hs, :a_v, 3)

# Couplings
@variables gu::Number gv::Number

###################################################
# QuantumCumulants.jl section
###################################################

# SLH triplets
G_u = SLH(1, gu*aun, 0) # input cavity
G_2lvls = SLH(1, √(γ)*σ(1, 2), 0) # 2-level system 
G_v = SLH(1, gv*avn, 0) # output cavity 

# I-O cascade
G_cas_1m = ▷(G_u, G_2lvls, G_v)

# Hamiltonian and Lindbladian
Hcas_1m = hamiltonian(G_cas_1m)
Lcas_1m = lindblad(G_cas_1m)[1]
Lcas_1md = adjoint(Lcas_1m)

# Time-dependent couplings
@register_symbolic gu_t(t)
@register_symbolic gv_t(t)
# Complex conjugates
@register_symbolic gu_c_t(t)
@register_symbolic gv_c_t(t)

gs_ls = [gu, gv]
gst_ls = [gu_t(t), gv_t(t)]
gsct_ls = [gu_c_t(t), gv_c_t(t)]

dict_gts = Dict([gs_ls; conj.(gs_ls)] .=> [gst_ls; gsct_ls]);

# Insert time-dependence
Hcas_1m_t = substitute(Hcas_1m, dict_gts);
Lcas_1m_t = substitute(Lcas_1m, dict_gts);
Lcas_1md_t = substitute(Lcas_1md, dict_gts);

0.5*1/(√(2π))*α1_/σ1_*exp(-(t-τ1_)^2/(2*σ1_^2))*exp(-1im*Δ1_*t)

# Pulses and composite mode
Ω1(t_) = α1_/(sqrt(2π*σ1_^2))*exp(-(t_-τ1_)^2/(2σ1_^2))*exp(-1im*Δ1_*t_)
Ω2(t_) = α2_/(sqrt(2π*σ2_^2))*exp(-(t_-τ2_)^2/(2σ2_^2))*exp(-1im*Δ2_*t_)
Ωall(t_) = Ω1(t_) + Ω2(t_)
αeff = sqrt(cumul_integrate(T, abs2.(Ωall.(T)))[end])
u(t_) = Ωall(t_) / αeff

gu_t_ = coupling_input(u, T)
gu_t(t) = gu_t_(t)
gv_t_ = coupling_output(u, T)
gv_t(t) = gv_t_(t)

gu_c_t(t) = conj(gu_t(t))
gv_c_t(t) = conj(gv_t(t))

# Second-order expansion for correlation-based mode extraction
order = 2
# Minimal operator set to capture second-order dynamics
ops_1m = [aun, aun'aun, aun*aun, σ(1, 2), σ(2, 2), avn, avn'avn]
# ops_all1 = unique([ops1[i]*ops1[j] for i = 1:length(ops1) for j = 1:length(ops1)]); # TODO: needed for corr?
# ops_all = union(ops1,ops_all1[1:1],ops_all1[3:14], ops_all1[16:16], ops_all1[18:26], ops_all1[28:29])

eqs_1m =
    meanfield(ops_1m, Hcas_1m_t, [Lcas_1m_t]; Jdagger = [Lcas_1md_t], order = order, iv = t);
complete!(eqs_1m);
length(eqs_1m)

# Classical-to-coherent relation
α_1m = αeff/(2*sqrt(γ_)) # may change if cross-terms are significant

u0_1m = zeros(ComplexF64, length(eqs_1m))
u0_1m[1] = α_1m
u0_1m[2] = abs2(α_1m)
u0_1m[3] = α_1m*α_1m

# Solve ODE system
sys_1m = mtkcompile(System(eqs_1m; name = :sys_1m))
ps = [γ]
p0 = [γ_]
u0_p_maps_1m = Dict([unknowns(sys_1m); γ] .=> [u0_1m; γ_])
prob_1m = ODEProblem(sys_1m, u0_p_maps_1m, (dt, Tend))
sol_1m = solve(prob_1m, Tsit5(), saveat = dt; abstol, reltol);

t_1m = sol_1m.t;
s22_1m = solution_values(sol_1m, σ(2, 2), eqs_1m);
nu_1m = real.(solution_values(sol_1m, aun'aun, eqs_1m))
nv_1m = real.(solution_values(sol_1m, avn'avn, eqs_1m))

common = (;
    xlims = (t_1m[1]-0.01, t_1m[end]),
    tickfontsize = plot_font,
    guidefontsize = plot_font,
    legendfontsize = plot_font,
)
p1 = plot(
    t_cl,
    s22_cl;
    color = :blue,
    label = L"\mathrm{classical}",
    ylabel = L"P",
    xticks = ([0, 10, 20], ["", "", ""]),
    legend = :topleft,
    common...,
)
plot!(p1, t_1m, real.(s22_1m); color = :red, ls = :dash, label = L"\mathrm{single~SUPER}")
p2 = plot(
    t_1m,
    nu_1m;
    color = :blue,
    label = L"\langle \hat{n}_{u} \rangle",
    ylabel = L"\langle\hat{n}_i\rangle",
    xticks = ([0, 10, 20], ["", "", ""]),
    ylims = (nv_1m[1]-0.5e3, nu_1m[1]+0.5e3),
    legend = :left,
    common...,
)
plot!(p2, t_1m, nv_1m; color = :red, label = L"\langle \hat{n}_{v} \rangle")
p3 = plot(
    t_1m,
    nu_1m .+ nv_1m .- nu_1m[1];
    color = :blue,
    label = L"\mathrm{mode}",
    ylabel = L"\langle\Delta\hat{n}_i\rangle",
    xlabel = L"\gamma t",
    yticks = ([-1, -0.5, 0], latexstring.([-1, -0.5, 0])),
    legend = :bottomleft,
    common...,
)
pl4 = plot(p1, p2, p3; layout = (3, 1), size = (800, 800))
display(pl4)
# savefig(pl4, "plot_IO_cascade.pdf")


###################################################
########### I-O SUPER Correlation #################
###################################################

# Same as above but with sparser saveat (va = 0.05) and only Input–System cascade

dt = 1e-5
Tend = 12*2
T = [dt:dt:Tend;]
ΔT = T[2] - T[1]

@register_symbolic gu_t(t)
@register_symbolic gv_t(t)
@register_symbolic gu_c_t(t)
@register_symbolic gv_c_t(t)

G_cas_1m0 = ▷(G_u, G_2lvls, G_v)

Hcas_1m0 = hamiltonian(G_cas_1m0)
Lcas_1m0 = lindblad(G_cas_1m0)[1]
Lcas_1m0d = adjoint(Lcas_1m0)

# Set v-coupling to zero (collect only input-system correlations)
gs0_ls = [gu, gv]
gs0t_ls = [gu_t(t), 0.0*gv_t(t)]
gsc0t_ls = [gu_c_t(t), 0.0*gv_c_t(t)]

dict_gts0 = Dict([gs0_ls; conj.(gs0_ls)] .=> [gs0t_ls; gsc0t_ls]);

Hcas_1m0_t = substitute(Hcas_1m0, dict_gts0)
Lcas_1m0_t = substitute(Lcas_1m0, dict_gts0)
Lcas_1m0d_t = substitute(Lcas_1m0d, dict_gts0)

gu_t_ = coupling_input(u, T)
gu_t(t) = gu_t_(t)
gv_t_ = coupling_output(u, T)
gv_t(t) = 0.0*gv_t_(t)

gu_c_t(t) = conj(gu_t(t))
gv_c_t(t) = conj(gv_t(t))

# order = 2 
# ops1 = [aun, aun', σ(1,2), σ(2,2),σ(1,2)', avn, avn']
# ops_all1 = unique([ops1[i]*ops1[j] for i = 1:length(ops1) for j = 1:length(ops1)]);
# ops_all = union(ops1,ops_all1[1:1],ops_all1[3:14], ops_all1[16:16], ops_all1[18:26], ops_all1[28:29])
# eqs_1m0 = meanfield(ops_all, Hcas_1m0_t, [Lcas_1m0_t]; Jdagger=[Lcas_1m0d_t], order=order, iv=t)

# complete!(eqs_1m0)
# length(eqs_1m0)

ops_1m0 = [aun, aun'aun, aun*aun, σ(1, 2), σ(2, 2)]
eqs_1m0 = meanfield(
    ops_1m0,
    Hcas_1m0_t,
    [Lcas_1m0_t];
    Jdagger = [Lcas_1m0d_t],
    order = order,
    iv = t,
);
complete!(eqs_1m0);
length(eqs_1m0)

u0_1m0 = zeros(ComplexF64, length(eqs_1m0))
u0_1m0[1] = α_1m
u0_1m0[2] = abs2(α_1m)
u0_1m0[3] = α_1m*α_1m

va = 0.05

# Solve ODE system
sys_1m0 = mtkcompile(System(eqs_1m0; name = :sys_1m0))
u0_p_maps_1m0 = Dict([unknowns(sys_1m0); γ] .=> [u0_1m0; γ_])
prob_1m0 = ODEProblem(sys_1m0, u0_p_maps_1m0, (dt, Tend))
sol_1m0 = solve(prob_1m0, Tsit5(), saveat = va; abstol, reltol);

###################################################
########### Correlation Part ######################
###################################################

# Lindblad vector for correlation functions
Lvec=[aun, σ(1, 2)]

# Discrete time grid
N = length(sol_1m0.t)

# g^(1)(t1,t2) container for each operator pair
Gmat_bar = zeros(ComplexF64, N, length(Lvec), length(Lvec), N)
t1 = sol_1m0.t;
t2 = sol_1m0.t;
τ0 = 0.0

# All operator correlation combinations
Corrmat = [
    [
        CorrelationFunction(Lvec[k]', Lvec[l], eqs_1m0; steady_state = false, iv0 = t)
        for l = 1:length(Lvec)
    ] for k = 1:length(Lvec)
]
# Minimal set containing full dynamics
Corrvec = [Corrmat[end][i] for i = 1:length(Lvec)]

xx = CorrelationFunction(σ(1, 2), σ(1, 2), eqs_1m0; steady_state = false, iv0 = t)
xx.eqs.states

###################################################
###### QC.jl extended correlation utilities #######
###################################################

# Structural equality of two lhs expressions
function same_lhs(a, b)
    operation(a) === operation(b) && begin
        aa = Tuple(arguments(a))
        bb = Tuple(arguments(b))
        length(aa) == length(bb) && all(isequal.(aa, bb))
    end
end

# Index mapping between Corrvec and Corrmat
function idx_tuple(Lvec, Corrvec, Corrmat)
    idxvec = zeros(Int8, length(Lvec), length(Lvec))
    for l = 1:length(Lvec)
        eqs1 = equations(Corrvec[l].de)
        for k = 1:length(Lvec)
            target_lhs = equations(Corrmat[k][l].de)[1].lhs
            idx = findfirst(eq -> same_lhs(eq.lhs, target_lhs), eqs1)
            idxvec[k, l] = idx
        end
    end
    idx_tuple = Tuple(ntuple(i -> c[i], length(Lvec)) for c in eachcol(idxvec))
    return idx_tuple
end

# Operators included
oper_tup = (Lvec[1], Lvec[2])
# Index tuple for selecting correlation components
idx_tup = idx_tuple(Lvec, Corrvec, Corrmat)

# Compute full Gmat_bar on the discrete time grid
function compute_Gmat(
    sol_1m0,
    eqs_1m0,
    Lvec,
    Corrvec;
    p0,
    ps,
    saveat,
    τ0 = 0.0,
    alg = OrdinaryDiffEq.Tsit5(),
    state_syms = oper_tup,
    extract_idxs = idx_tup,
)

    N = length(sol_1m0.t)
    Gmat_bar = zeros(ComplexF64, N, length(Lvec), length(Lvec), N)

    padleftN(v::AbstractVector{<:Number}, N::Int) = (
        length(v) < N ? vcat(zeros(ComplexF64, N - length(v)), ComplexF64.(v)) :
        ComplexF64.(v)
    )

    for i1 = 1:(N-1)
        τend = sol_1m0.t[end] - sol_1m0.t[i1]
        z = i1

        for l = 1:length(Lvec)
            corr = Corrvec[l]
            complete(corr.de)

            iv = ModelingToolkit.get_iv(corr.de)
            de_shift = substitute(corr.de, Dict(t => iv + sol_1m0.t[z]))

            bc_var = complete(corr.de)[end].lhs
            rhs_val = sol_1m0[state_syms[l]][z]
            corrn_fin = substitute(de_shift, Dict(bc_var => rhs_val))

            @named csys = System(corrn_fin)

            u0_c = correlation_u0(corr, sol_1m0.u[z])
            p0_c = correlation_p0(corr, sol_1m0.u[z], ps .=> p0)
            dict = merge(Dict(u0_c), Dict(p0_c))

            prob_c = ODEProblem(csys, dict, (τ0, τend))
            sol_c = solve(prob_c, alg; saveat = saveat)

            idx_v = extract_idxs[l]
            u_flat = [ComplexF64.(getindex.(sol_c.u, idx_v[i])) for i = 1:length(Lvec)]

            if i1 == 1
                for k = 1:length(Lvec)
                    Gmat_bar[i1, k, l, :] .= u_flat[k]
                end
            else
                for k = 1:length(Lvec)
                    Gmat_bar[i1, k, l, :] .= padleftN(u_flat[k], N)
                end
            end
        end
    end

    return Gmat_bar
end


# Optionally compute and save Gmat_bar
#Gmat_bar = compute_Gmat(
#    sol_1m0, eqs_1m0, Lvec, Corrvec;
#    p0 = p0,
#    ps = ps,
#    saveat = va,
#    τ0 = 0.0,
#    alg = Tsit5(),
#    state_syms = oper_tup,
#    extract_idxs = idx_tup
#)

# Load precomputed correlation data
#JLD2.save_object("Corr_Single.jld2",[Gmat_bar,t1])
Gmat_bar = JLD2.load_object("Corr_Single.jld2")[1];
t1 = JLD2.load_object("Corr_Single.jld2")[2];
N = length(t1)
# Build g^(1)(t1,t2)
gcorr1 = zeros(ComplexF64, N, N)
gcorr = zeros(ComplexF64, N, N)

for i1 = 1:N
    for j1 = 1:N
        gcorr1[j1, i1] =
            gu_c_t(t1[j1])*gu_t(t1[i1])*Gmat_bar[i1, 1, 1, j1] +
            gu_c_t(t1[j1])*sqrt(γ_)*Gmat_bar[i1, 1, 2, j1] +
            gu_t(t1[i1])*sqrt(γ_)*Gmat_bar[i1, 2, 1, j1] +
            γ_*Gmat_bar[i1, 2, 2, j1]
    end
end
gcorr1
Λ = zeros(ComplexF64, N, N)
for i = 1:N
    Λ[i, i] = diag(real.(gcorr1))[i]
end

# Symmetrized final g^(1)
gcorr = gcorr1 + tril(gcorr1, 1)' - Λ


###################################################
########### Plots: I-O SUPER Corr. ################
###################################################

###################################################
# Heat plot
###################################################

# Normalization for visualization
g_abs = abs.(gcorr)
g_max = maximum(g_abs)
g_vis = g_abs ./ (g_max > 0 ? g_max : 1)

tick_vals = collect(0.0:0.2:1.0)
tick_labs = latexstring.(string.(round.(tick_vals; digits = 1)))
pl6 = heatmap(
    t1,
    t1,
    g_vis;
    color = :thermal,
    clims = (0, 1),
    xlabel = L"\gamma t_2",
    ylabel = L"\gamma t_1",
    colorbar_title = L"|g^{(1)}|\ \mathrm{(norm)}",
    colorbar_ticks = (tick_vals, tick_labs),
    aspect_ratio = :equal,
    tickfontsize = 16,
    guidefontsize = 18,
    size = (600, 500),
)
display(pl6)
# savefig(pl6, "heat_IO.pdf")


#####################################################################
# Eigenvalues and Eigenmodes
#####################################################################

F = eigen(gcorr);

ΔT = t1[2]-t1[1]
ns_avg = F.values
modess = F.vectors
vs_mode = (modess[:, end]) / sqrt(ΔT)
ns_avgt = ns_avg*(t1[2] - t1[1])
ns = ns_avgt[end]
ns_avgt

###################################################
########### Plot I-O SUPER Cascade ################
###################################################

# Dominant mode magnitude vs composite input
pl7 = plot(
    t1,
    abs.(u.(t1));
    color = :blue,
    label = L"v(t) = u(t)",
    ylabel = L"|v(t)|",
    xlims = (t1[1]-0.01, t1[end]),
    xticks = ([0, 10, 20], ["", "", ""]),
    legend = :topleft,
    tickfontsize = plot_font,
    guidefontsize = plot_font,
    legendfontsize = plot_font,
    size = (800, 200),
)
plot!(pl7, t1, abs.(vs_mode); color = :red, ls = :dash, label = L"\mathrm{corr.}")
display(pl7)
# savefig(pl7, "mode_IO.pdf")


###################################################
## Loss mode: I-O SUPER with corr-based approach ##
###################################################
dt = 1e-5
Tend = 12*2
T = [dt:dt:Tend;]
ΔT = T[2] - T[1]

# Hilbert space
hu = FockSpace(:u) # virtual input cavity 
hv = FockSpace(:v) # virtual output cavity 
hs = tensor(hu, hs1, hv)

# Operators
aun = Destroy(hs, :a_u, 1)
σ(i, j) = Transition(hs, :ss, i, j)
avn = Destroy(hs, :a_v, 3)

# Couplings
gu, gv = cnumbers("g_{u} g_{v}")


###################################################
# QuantumCumulants.jl section
###################################################

# SLH triplets
G_u = SLH(1, gu*aun, 0) # input cavity
G_2lvls = SLH(1, √(γ)*σ(1, 2), 0) # 2-level system 
G_v = SLH(1, gv*avn, 0) # output cavity 

# I-O cascade
G_cas_1m = ▷(G_u, G_2lvls, G_v)

# Hamiltonian and Lindbladian
Hcas_1m = hamiltonian(G_cas_1m)
Lcas_1m = lindblad(G_cas_1m)[1]
Lcas_1md = adjoint(Lcas_1m)

# Time-dependent couplings
@register_symbolic gu_t(t)
@register_symbolic gv_t(t)
# Complex conjugates
@register_symbolic gu_c_t(t)
@register_symbolic gv_c_t(t)

gs_ls = [gu, gv]
gst_ls = [gu_t(t), gv_t(t)]
gsct_ls = [gu_c_t(t), gv_c_t(t)]

dict_gts = Dict([gs_ls; conj.(gs_ls)] .=> [gst_ls; gsct_ls]);

# Insert time-dependence
Hcas_1m_t = substitute(Hcas_1m, dict_gts);
Lcas_1m_t = substitute(Lcas_1m, dict_gts);
Lcas_1md_t = substitute(Lcas_1md, dict_gts);

# Pulses and composite mode
αeff = 0.0

Ω1(t_) = α1_/(sqrt(2π*σ1_^2)) * exp(-(t_-τ1_)^2/(2σ1_^2)) * exp(-1im*Δ1_*t_)

Ω2(t_) = α2_/(sqrt(2π*σ2_^2)) * exp(-(t_-τ2_)^2/(2σ2_^2)) * exp(-1im*Δ2_*t_)

Ωall(t_) = Ω1(t_) + Ω2(t_)

αeff = sqrt(cumul_integrate(T, abs2.(Ωall.(T)))[end])

u(t_) = Ωall(t_) / αeff

gu_t_ = coupling_input(u, T)
gu_t(t) = gu_t_(t)
gv_t_ = coupling_output(u, T)
gv_t(t) = gv_t_(t)


gu_c_t(t) = conj(gu_t(t))
gv_c_t(t) = conj(gv_t(t))


order = 2
ops1 = [aun, aun', σ(1, 2), σ(2, 2), σ(1, 2)', avn, avn']
ops_all1 = unique([ops1[i]*ops1[j] for i = 1:length(ops1) for j = 1:length(ops1)]);
ops_all = union(
    ops1,
    ops_all1[1:1],
    ops_all1[3:14],
    ops_all1[16:16],
    ops_all1[18:26],
    ops_all1[28:29],
)
eqs_1m = meanfield(
    ops_all,
    Hcas_1m_t,
    [Lcas_1m_t];
    Jdagger = [Lcas_1md_t],
    order = order,
    iv = t,
)

complete!(eqs_1m)
length(eqs_1m)

# Classical-to-coherent relation
α_1m = αeff/(2*sqrt(γ_)) # may change if cross-terms are significant

function ψ0_us1v(α, eqs)
    bu = FockBasis(rd(abs2(α) + 20abs(α)), rd(abs2(α) - 20abs(α))) # include ±20α photon fluctuations
    bs1 = NLevelBasis(2)
    bv = FockBasis(1)
    b = tensor([bu, bs1, bv]...)
    ψu = coherentstate(bu, α)
    ψs1 = nlevelstate(bs1, 1)
    ψv = fockstate(bv, 0)
    ψ0 = LazyKet(b, (ψu, ψs1, ψv))
    return initial_values(eqs, ψ0)
end

# Initial values
u0_1m = ψ0_us1v(α_1m, eqs_1m)

va = 0.05

# Solve ODE system
sys_1m = mtkcompile(System(eqs_1m; name = :sys_1m))
ps = [γ]
p0 = [γ_]
u0maps = Dict(unknowns(sys_1m) .=> u0_1m)
pmaps = parameter_map(sys_1m, ps, p0)
probs = ODEProblem(sys_1m, merge(u0maps, pmaps), (dt, Tend))
sol_1m0v = solve(probs, Tsit5(), saveat = va; abstol, reltol);


###################################################
########### Correlation Part ######################
###################################################

# Lindblad vector (includes output)
Lvec=[aun, σ(1, 2), avn]

N = length(sol_1m0v.t)
Gmat_bar = zeros(ComplexF64, N, length(Lvec), length(Lvec), N)
t1 = sol_1m0v.t;
t2 = sol_1m0v.t;
τ0 = 0.0

Corrmat = [
    [
        CorrelationFunction(Lvec[k]', Lvec[l], eqs_1m; steady_state = false) for
        l = 1:length(Lvec)
    ] for k = 1:length(Lvec)
]
Corrvec = [Corrmat[end][i] for i = 1:length(Lvec)]

oper_tup = (Lvec[1], Lvec[2], Lvec[3])
idx_tup = idx_tuple(Lvec, Corrvec, Corrmat)

# Optionally compute and save
#Gmat_bar = compute_Gmat( 
#    sol_1m0v, eqs_1m, Lvec, Corrvec;
#    p0 = p0,
#    ps = ps,
#    saveat = va,
#    τ0 = 0.0,
#    alg = Tsit5(),
#    state_syms = oper_tup,
#    extract_idxs = idx_tup
#)

# Load precomputed data
#JLD2.save_object("Corr_v.jld2",[Gmat_bar,t1])
Gmat_bar = JLD2.load_object("Corr_v.jld2")[1];
t1 = JLD2.load_object("Corr_v.jld2")[2];
N = length(t1)
gcorr1 = zeros(ComplexF64, N, N)
gcorr = zeros(ComplexF64, N, N)

gvec(t) = [gu_t(t), sqrt(γ_), gv_t(t)]
gvecc(t) = [gu_c_t(t), sqrt(γ_), gv_c_t(t)]

for i1 = 1:N
    for j1 = 1:N
        gcorr1[i1, j1] = sum(
            gvecc(t1[i1])[k]*gvec(t1[j1])[l]*(Gmat_bar[j1, k, l, i1]) for k = 1:3 for
            l = 1:3
        )
    end
end
gcorr1
Λ = zeros(ComplexF64, N, N)
for i = 1:N
    Λ[i, i] = diag(real.(gcorr1))[i]
end

# Final g^(1)(t1,t2)
gcorr = gcorr1 + tril(gcorr1, 1)' - Λ



###################################################
########### Plot I-O SUPER Corr. loss #############
###################################################

####################################################################
# Heat plot
####################################################################

g_abs = abs.(gcorr)
g_max = maximum(g_abs)
g_vis = g_abs ./ (g_max > 0 ? g_max : 1)

tick_vals = collect(0.0:0.2:1.0)
tick_labs = latexstring.(string.(round.(tick_vals; digits = 1)))
pl8 = heatmap(
    t1,
    t1,
    g_vis;
    color = :thermal,
    clims = (0, 1),
    xlabel = L"\gamma t_2",
    ylabel = L"\gamma t_1",
    colorbar_title = L"|g^{(1)}|\ \mathrm{(norm)}",
    colorbar_ticks = (tick_vals, tick_labs),
    aspect_ratio = :equal,
    tickfontsize = 16,
    guidefontsize = 18,
    size = (600, 500),
)
display(pl8)
# savefig(pl8, "heat_IO_loss.pdf")


#####################################################################
# Eigenvalues and eigenmodes
#####################################################################

F = eigen(gcorr);

ΔT = t1[2]-t1[1]
ns_avg = F.values
modess = F.vectors
vs_mode = (modess[:, end]) / sqrt(ΔT)
ns_avgt = ns_avg*(t1[2] - t1[1])
ns = ns_avgt[end]
ns_avgt

###################################################
########### Plot I-O SUPER Cascade mode ###########
###################################################

# Dominant mode magnitude (loss-inclusive)
pl9 = plot(
    t1,
    abs.(u.(t1));
    color = :blue,
    label = L"v(t) = u(t)",
    ylabel = L"|v(t)|",
    xlims = (t1[1]-0.01, t1[end]),
    xticks = ([0, 10, 20], ["", "", ""]),
    legend = :topleft,
    tickfontsize = plot_font,
    guidefontsize = plot_font,
    legendfontsize = plot_font,
    size = (800, 200),
)
plot!(pl9, t1, abs.(vs_mode); label = L"\mathrm{loss}")
display(pl9)
# savefig(pl9, "mode_IO_loss.pdf")
