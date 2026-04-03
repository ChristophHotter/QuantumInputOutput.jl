###########################
### interaction picture ###
###########################

"""
    interaction_picture_M(A, T; alg = Tsit5(), kwargs...)

Solve the Heisenberg coefficient-matrix equation for an interaction picture,
`dM/dt = A(t) * M(t)` with `M(0) = I`, and return the ODE solution.

The function `A(t)` must return an `N × N` complex matrix describing the linear
mode coupling in the interaction picture. The returned solution supports
interpolation `M(t)` for continuous `t`.
"""
function interaction_picture_M(A::Function, T; alg = Tsit5(), kwargs...)
    T0 = T[1]
    Tend = T[end]
    n = size(A(T0), 1)
    M0 = Matrix{ComplexF64}(I, n, n)
    f_M(u, p, t) = A(t) * u
    prob_M = ODEProblem(f_M, M0, (T0, Tend))
    sol = solve(prob_M, alg; saveat = T, kwargs...)
    return t -> sol(t)
end

"""
    interaction_picture_M_2modes_equal(u, T)

Analytic interaction-picture coefficient matrix for two modes when `u(t) = v(t)`.
The matrix is

```math
M(t) = \\begin{bmatrix}
\\cos \\theta(t) & -\\sin \\theta(t) \\\\
\\sin \\theta(t) & \\cos \\theta(t)
\\end{bmatrix},
```

where

```math
\\sin^2 \\theta(t) = \\int_0^t |u(t')|^2\\,dt'.
```
"""
function interaction_picture_M_2modes_equal(u, T)
    u_vals = u isa Function ? u.(T) : u
    sin2θ = cumul_integrate(T, abs2.(u_vals))
    sin2θ = clamp.(real.(sin2θ), 0.0, 1.0)
    θ = asin.(sqrt.(sin2θ))
    cθ = cos.(θ)
    sθ = sin.(θ)

    M11 = LinearInterpolation(cθ, T; extrapolation = ExtrapolationType.Extension)
    M12 = LinearInterpolation(-sθ, T; extrapolation = ExtrapolationType.Extension)
    M21 = LinearInterpolation(sθ, T; extrapolation = ExtrapolationType.Extension)
    M22 = LinearInterpolation(cθ, T; extrapolation = ExtrapolationType.Extension)

    return t -> [M11(t) M12(t); M21(t) M22(t)]
end

_as_time_function(x) = x isa Function ? x : (_ -> x)

"""
    interaction_picture_A(gs)
    interaction_picture_A(g1, g2, gs...)

Coefficient matrix `A(t)` for `N` interacting modes with couplings `gs = [g_1, …, g_N]`.
All couplings may be time-dependent functions or constants.

The anti-Hermitian matrix has entries

```math
A_{ij}(t) = \\tfrac{1}{2}\\,g_i(t)\\,g_j^*(t), \\quad i < j,
\\qquad A_{ji} = -A_{ij}^*,
\\qquad A_{ii} = 0.
```

For two modes this gives

```math
A(t) = \\frac{1}{2}
\\begin{bmatrix}
0 & g_1(t)\\, g_2^*(t) \\\\
-g_1^*(t)\\, g_2(t) & 0
\\end{bmatrix},
```

and for three modes

```math
A(t) = \\frac{1}{2}
\\begin{bmatrix}
0 & g_1(t)\\, g_2^*(t) & g_1(t)\\, g_3^*(t) \\\\
-g_1^*(t)\\, g_2(t) & 0 & g_2(t)\\, g_3^*(t) \\\\
-g_1^*(t)\\, g_3(t) & -g_2^*(t)\\, g_3(t) & 0
\\end{bmatrix}.
```
"""
function interaction_picture_A(gs)
    gfs = _as_time_function.(gs)
    n = length(gfs)
    function A(t)
        g = [gf(t) for gf in gfs]
        M = zeros(ComplexF64, n, n)
        for i in 1:n, j in (i+1):n
            M[i, j] = g[i] * conj(g[j])
            M[j, i] = -conj(M[i, j])
        end
        return 0.5 * M
    end
    return A
end

interaction_picture_A(g1, g2, gs...) = interaction_picture_A([g1, g2, gs...])
