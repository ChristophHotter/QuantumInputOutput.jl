# # Traveling Cat State from a Kerr-Nonlinear Parametric Oscillator

# In this example, we reproduce the Wigner function shown in Fig. 3(b) of
# [A. Kiilerich and K. Molmer, Phys. Rev. A 102, 023717 (2020)](https://doi.org/10.1103/PhysRevA.102.023717).
#We drive a Kerr-nonlinear parametric
# oscillator (KPO), extract the dominant output mode from the cavity-field
# autocorrelation function, and then capture this mode in a virtual output
# cavity. The Wigner function at the $t = 31 K^{-1}$ is plotted.

# We start by loading the packages and defining the symbolic operators and
# parameters for the driven KPO.

using QuantumInputOutput
using SecondQuantizedAlgebra
using QuantumOptics
using Plots
using LinearAlgebra
using DataInterpolations

#

## symbolic Hilbert space
hc = FockSpace(:c)

## symbolic operator
a = Destroy(hc, :a, 1)

## symbolic parameters
@rnumbers γ K Δ
p = cnumber("p")
nothing # hide

# The KPO Hamiltonian is
# $H_s(t) = \frac{p(t)}{2}(a^{\dagger 2} + a^2) - \frac{K}{2} a^{\dagger 2} a^2 + \Delta a^\dagger a$
# and the cavity emits into the waveguide with jump operator $\sqrt{\gamma} a$.

H_s = p / 2 * (a'^2 + a^2) - K / 2 * (a'^2) * (a^2) + Δ * a' * a
G_s = SLH(1, √(γ) * a, H_s)

H = get_hamiltonian(G_s)
L = get_lindblad(G_s)[1]
nothing # hide

# Next, we define the numerical parameters from Sec. II.D. The classical pump
# is obtained from a fourth-order low-pass filter driven by
# $p_{in}(t) = K A_p e^{-\gamma t}$.

γ_ = 1.0
K_ = 5γ_
Δ_ = 0.0
B_ = 2.5γ_
A_p = 4.45

T_end = 14 / γ_
T = [0:0.004:1;]*T_end
ΔT = T[2] - T[1]
nothing # hide

function lowpass_filter(signal, T, B)
    filtered = zeros(ComplexF64, length(signal))
    for i = 2:length(T)
        dt = T[i] - T[i - 1]
        α = exp(-B * dt)
        filtered[i] = α * filtered[i - 1] + (1 - α) * signal[i - 1]
    end
    return filtered
end

pin = ComplexF64[K_ * A_p * exp(-γ_ * t) for t in T]
pump = copy(pin)
for _ = 1:4
    pump = lowpass_filter(pump, T, B_)
end
p_t_ = LinearInterpolation(pump, T)
p_t(t) = p_t_(t)
nothing # hide

# We solve the master equation for the driven KPO and use the first-order
# correlation matrix $g^{(1)}(t_1, t_2)$ to determine the dominant output mode.

n_cut = 10
bc = FockBasis(n_cut)
b1 = bc

dict_p = Dict([γ, K, Δ] .=> [γ_, K_, Δ_])
dict_p_t = Dict(p => p_t)

H_QO = translate_qo(H, b1; parameter = dict_p, time_parameter = dict_p_t)
L_QO = translate_qo(L, b1; parameter = dict_p, time_parameter = dict_p_t)

function input_output_1(t, ρ)
    Ht = H_QO(t)
    J = [L_QO(t)]
    return Ht, J, dagger.(J)
end

ψ0 = fockstate(bc, 0)
t_1, ρt_1 = timeevolution.master_dynamic(T, ψ0, input_output_1)
nothing # hide

a_qo = destroy(bc)
Ls(t) = √(γ_) * a_qo
g1_m = two_time_corr_matrix(T, ρt_1, input_output_1, Ls)
# g1_m = two_time_corr_matrix(T, ρt_1, input_output_1, √(γ_) * a_qo) # TODO: add method!

F = eigen(g1_m)
n_avg = real.(F.values) * ΔT
v_mode = F.vectors[:, end] / √(ΔT)
nothing # hide

# After identifying the dominant output mode, we add a virtual output cavity
# that captures this mode explicitly.

hv = FockSpace(:v)
h = hc ⊗ hv

a2 = Destroy(h, :a, 1)
av = Destroy(h, :a_v, 2)
gv = cnumber("g_v")

H_s2 = p / 2 * (a2'^2 + a2^2) - K / 2 * (a2'^2) * (a2^2) + Δ * a2' * a2
G_s2 = SLH(1, √(γ) * a2, H_s2)
G_v = SLH(1, gv' * av, 0)
G = G_s2 ▷ G_v

H_2 = get_hamiltonian(G)
L_2 = get_lindblad(G)[1]

gv_t = v_to_gv(v_mode, T)
dict_p_t_2 = Dict(p => p_t, gv => gv_t)

bv = FockBasis(n_cut)
b2 = bc ⊗ bv

H_QO_2 = translate_qo(H_2, b2; parameter = dict_p, time_parameter = dict_p_t_2)
L_QO_2 = translate_qo(L_2, b2; parameter = dict_p, time_parameter = dict_p_t_2)

function input_output_2(t, ρ)
    Ht = H_QO_2(t)
    J = [L_QO_2(t)]
    return Ht, J, dagger.(J)
end

ψ0_2 = fockstate(bc, 0) ⊗ fockstate(bv, 0)
t_2, ρt_2 = timeevolution.master_dynamic(T, ψ0_2, input_output_2)
nothing # hide

# We monitor the KPO excitation, the dominant output-mode envelope, 
# and the fidelity of the captured state with the target cat state.

a2_qo = destroy(bc) ⊗ one(bv)
n_c_t = real.(expect(a2_qo' * a2_qo, ρt_2))
mode_shape = abs.(v_mode) 
nothing # hide

# TODO: look at this!!

# The output-mode eigenvector is only defined up to a global phase. To match
# the convention used in Fig. 3(b), we rotate the captured state such that the
# phase of $\langle a_v^2 \rangle$ agrees with the cat-state amplitude quoted
# in the paper, $\beta = 2 e^{-0.31 i \pi}$.

t_target = 31 / K_
target_index = argmin(abs.(T .- t_target))

ρv_target = ptrace(ρt_2[target_index], 1)
a_v = destroy(bv)

β_target = 2 * exp(-0.31im * π)
m2_target = expect(a_v^2, ρv_target)
phase_shift = (angle(m2_target) - 2 * angle(β_target)) / 2
U_phase = diagonaloperator(bv, exp.(-1im .* phase_shift .* (0:n_cut)))
ρv_plot = U_phase * ρv_target * dagger(U_phase)

cat_state = normalize(coherentstate(bv, β_target) + coherentstate(bv, -β_target))
cat_proj = projector(cat_state)
ρv_t_rot = [U_phase * ptrace(ρ, 1) * dagger(U_phase) for ρ in ρt_2]
fidelity_t = [real(expect(cat_proj, ρv)) for ρv in ρv_t_rot]
nothing # hide

#

p1 = plot(
    T,
    fidelity_t;
    color = :red,
    lw = 2,
    label = "cat fidelity",
)
plot!(
    p1,
    T,
    mode_shape;
    color = :grey,
    fillrange = 0,
    fillalpha = 0.25,
    label = "mode shape",
)
plot!(
    p1;
    xlabel = "time (1/γ)",
    ylabel = "fidelity / mode shape",
    xlims = (0, T_end),
    ylims = (0, 1),
    legend = :topright,
    size = (560, 320),
)

p1r = twinx(p1)
plot!(
    p1r,
    T,
    n_c_t;
    color = :blue,
    lw = 2,
    ls = :dot,
    label = "KPO excitation",
)
plot!(p1r; ylabel = "KPO excitation")
p1

# Finally, we plot only the Wigner function of the captured traveling mode at
# $t = 31 K^{-1}$.

x = [-5:0.05:5;]
W = wigner(ρv_plot, x, x)
W_max = maximum(abs, W)

p = heatmap(
    x/√(2),
    x/√(2),
    real.(W');
    c = cgrad(:RdBu_11, rev = true),
    clims = (-W_max, W_max),
    xlabel = "x",
    ylabel = "p",
    title = "t = 31/K",
    xlims = (-3.5, 3.5),
    ylims = (-3.5, 3.5),
    aspect_ratio = 1,
    size = (460, 380),
)
p

# ## Package versions

# These results were obtained using the following versions:

using InteractiveUtils
versioninfo()

using Pkg
Pkg.status(
    [
        "QuantumInputOutput",
        "SecondQuantizedAlgebra",
        "QuantumOptics",
        "Plots",
        "DataInterpolations",
    ],
    mode = PKGMODE_MANIFEST,
)
