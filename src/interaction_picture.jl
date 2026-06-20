###########################
### interaction picture ###
###########################

using StaticArrays: SMatrix

_as_time_function(x::Number) = _ -> x
_as_time_function(x) = x  # anything callable passes through

"""
    coupling_matrix(gs::Tuple)

Build the antisymmetric coupling coefficient matrix `A(t)` from a tuple of
coupling functions/constants `gs = (g_1, ..., g_N)`. Returns a closure `t -> A(t)`.

```math
A_{ij}(t) = \\frac{1}{2} \\begin{cases}
0 & i = j \\\\
g_i(t)\\, g_j^*(t) & i < j \\\\
-g_i^*(t)\\, g_j(t) & i > j
\\end{cases}
```

All couplings may be time-dependent or constant.
"""
function coupling_matrix(gs::NTuple{N}) where {N}
    gfs = map(_as_time_function, gs)
    function A(t)
        gvals = ntuple(i -> gfs[i](t), Val(N))
        SMatrix{N,N,ComplexF64}(ntuple(Val(N * N)) do k
            i, j = divrem(k - 1, N) .+ (1, 1)
            if i == j
                zero(ComplexF64)
            elseif j < i
                0.5 * gvals[j] * conj(gvals[i])
            else  # j > i
                -conj(0.5 * gvals[i] * conj(gvals[j]))
            end
        end)
    end
    return A
end

coupling_matrix(g1, g2, gs...) = coupling_matrix((g1, g2, gs...))

"""
    solve_mode_evolution(A::Function, T; alg=Tsit5(), kwargs...)

Solve the interaction-picture coefficient-matrix ODE `dM/dt = A(t) M(t)` with `M(0) = I`.
Returns the ODE solution directly (callable as `sol(t)`).

All kwargs are passed on to the ODE solver.
"""
function solve_mode_evolution(A::Function, T; alg = Tsit5(), kwargs...)
    T0 = T[1]
    Tend = T[end]
    n = size(A(T0), 1)
    M0 = Matrix{ComplexF64}(I, n, n)
    function f_M!(du, u, p, t)
        mul!(du, A(t), u)
    end
    prob = ODEProblem(f_M!, M0, (T0, Tend))
    sol = solve(prob, alg; saveat = T, kwargs...)
    return sol
end

"""
    solve_mode_evolution_symmetric(u, T)

Analytic interaction-picture coefficient matrix for two modes when `u(t) = v(t)`.
Returns a callable `t -> M(t)` where

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
function solve_mode_evolution_symmetric(u, T)
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

    return t -> SMatrix{2,2}(M11(t), M21(t), M12(t), M22(t))  # column-major
end
