###########################
### interaction picture ###
###########################

"""
    interaction_picture_M(A, T; alg = Tsit5(), kwargs...)

Solve the Heisenberg coefficient-matrix equation for an interaction picture,
`dM/dt = A(t) * M(t)` with `M(T[1]) = I`, and return the ODE solution.

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

_as_time_function(x) = x isa Function ? x : (_ -> x)

"""
    interaction_picture_A_2modes(g1, g2)

Coefficient matrix `A(t)` for two virtual modes `(1, 2)`, where

`d/dt [a1; a2] = A(t) * [a1; a2]` and
`A(t) = 1/2 * [0, g1(t) * conj(g2(t)); -conj(g1(t)) * g2(t), 0]`.

All couplings may be time-dependent or constant.
"""
function interaction_picture_A_2modes(g1, g2)
    g1f = _as_time_function(g1)
    g2f = _as_time_function(g2)
    A(t_) = 0.5 * [0 g1f(t_) * conj(g2f(t_));
        -conj(g1f(t_)) * g2f(t_) 0]
    return A
end
interaction_picture_A_2modes(g_ls) = interaction_picture_A_2modes(g_ls...)

"""
    interaction_picture_A_3modes(g1, g2, g3)

Coefficient matrix `A(t)` for three interacting modes ordered as `(1, 2, 3)`,
with couplings `g1(t)`, `g2(t)` and `g3(t)` corresponding to Eq. (25) of
Christiansen et al. (PRA 107, 013706).

All couplings may be time-dependent or constant.
"""
function interaction_picture_A_3modes(g1, g2, g3)
    g1f = _as_time_function(g1)
    g2f = _as_time_function(g2)
    g3f = _as_time_function(g3)
    A(t_) = begin
        0.5 * [0 conj(g2f(t_)) * g1f(t_) g1f(t_) * conj(g3f(t_));
            -g2f(t_) * conj(g1f(t_)) 0 g2f(t_) * conj(g3f(t_));
            -conj(g1f(t_)) * g3f(t_) -conj(g2f(t_)) * g3f(t_) 0]
    end
    return A
end

"""
    interaction_picture_A_4modes(g1, g2, g3, g4)

Coefficient matrix `A(t)` for four interacting modes ordered as `(1, 2, 3, 4)`.
All couplings may be time-dependent or constant.
"""
function interaction_picture_A_4modes(g1, g2, g3, g4)
    g1f = _as_time_function(g1)
    g2f = _as_time_function(g2)
    g3f = _as_time_function(g3)
    g4f = _as_time_function(g4)
    A(t_) = 0.5 * [0 -g1f(t_) * conj(g2f(t_)) g1f(t_) * conj(g3f(t_)) g1f(t_) * conj(g4f(t_));
        conj(g1f(t_)) * g2f(t_) 0 g2f(t_) * conj(g3f(t_)) g2f(t_) * conj(g4f(t_));
        -conj(g1f(t_)) * g3f(t_) -conj(g2f(t_)) * g3f(t_) 0 g3f(t_) * conj(g4f(t_));
        -conj(g1f(t_)) * g4f(t_) -conj(g2f(t_)) * g4f(t_) -conj(g3f(t_)) * g4f(t_) 0]
    return A
end
