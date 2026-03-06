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
    interaction_picture_A_2modes(g1, g2)

Coefficient matrix `A(t)` for two virtual modes `(1, 2)`,

```math
A(t) = \\frac{1}{2}
\\begin{bmatrix}
0 & g_1(t) g_2^*(t) \\\\
-g_1^*(t) g_2(t) & 0
\\end{bmatrix}
```

All couplings may be time-dependent or constant.
"""
function interaction_picture_A_2modes(g1, g2)
    g1f = _as_time_function(g1)
    g2f = _as_time_function(g2)
    A(t) = 0.5 * [
        0 g1f(t) * conj(g2f(t));
        -conj(g1f(t)) * g2f(t) 0
    ]
    return A
end
interaction_picture_A_2modes(g_ls) = interaction_picture_A_2modes(g_ls...)

"""
    interaction_picture_A_3modes(g1, g2, g3)

Coefficient matrix `A(t)` for three interacting modes ordered as `(1, 2, 3)`

```math
A(t) = \\frac{1}{2}
\\begin{bmatrix}
0 & g_1(t) g_2^*(t) & g_1(t) g_3^*(t) \\\\
-g_1^*(t) g_2(t) & 0 & g_2(t) g_3^*(t) \\\\
-g_1^*(t) g_3(t) & -g_2^*(t) g_3(t) & 0
\\end{bmatrix}
```

All couplings may be time-dependent or constant.
"""
function interaction_picture_A_3modes(g1, g2, g3)
    g1f = _as_time_function(g1)
    g2f = _as_time_function(g2)
    g3f = _as_time_function(g3)
    A(t) =
        0.5 * [
            0 conj(g2f(t)) * g1f(t) g1f(t) * conj(g3f(t));
            -g2f(t) * conj(g1f(t)) 0 g2f(t) * conj(g3f(t));
            -conj(g1f(t)) * g3f(t) -conj(g2f(t)) * g3f(t) 0
        ]
    return A
end

"""
    interaction_picture_A_4modes(g1, g2, g3, g4)

Coefficient matrix `A(t)` for four interacting modes ordered as `(1, 2, 3, 4)`

```math
A(t) = \\frac{1}{2}
\\begin{bmatrix}
0 & g_1(t) g_2^*(t) & g_1(t) g_3^*(t) & g_1(t) g_4^*(t) \\\\
-g_1^*(t) g_2(t) & 0 & g_2(t) g_3^*(t) & g_2(t) g_4^*(t) \\\\
-g_1^*(t) g_3(t) & -g_2^*(t) g_3(t) & 0 & g_3(t) g_4^*(t) \\\\
-g_1^*(t) g_4(t) & -g_2^*(t) g_4(t) & -g_3^*(t) g_4(t) & 0
\\end{bmatrix}
```

All couplings may be time-dependent or constant.
"""
function interaction_picture_A_4modes(g1, g2, g3, g4)
    g1f = _as_time_function(g1)
    g2f = _as_time_function(g2)
    g3f = _as_time_function(g3)
    g4f = _as_time_function(g4)
    A(t) =
        0.5 * [
            0 g1f(t) * conj(g2f(t)) g1f(t) * conj(g3f(t)) g1f(t) * conj(g4f(t));
            -conj(g1f(t)) * g2f(t) 0 g2f(t) * conj(g3f(t)) g2f(t) * conj(g4f(t));
            -conj(g1f(t)) * g3f(t) -conj(g2f(t)) * g3f(t) 0 g3f(t) * conj(g4f(t));
            -conj(g1f(t)) * g4f(t) -conj(g2f(t)) * g4f(t) -conj(g3f(t)) * g4f(t) 0
        ]
    return A
end
