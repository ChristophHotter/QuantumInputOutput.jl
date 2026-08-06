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
-g_j^*(t)\\, g_i(t) & i > j
\\end{cases}
```

so that `A(t)` is anti-Hermitian. All couplings may be time-dependent or constant.
"""
function coupling_matrix(gs::Tuple{Vararg{Any,N}}) where {N}
    gfs = map(_as_time_function, gs)
    A(t) = _coupling_matrix(gfs, t)
    return A
end

# Unrolled at compile time so that tuples of mixed element type (a constant next to
# an interpolant, say) stay allocation-free: indexing such a tuple in a loop is only
# type stable if constant propagation reaches the index, and it does not.
@generated function _coupling_matrix(gfs::Tuple{Vararg{Any,N}}, t) where {N}
    g = [Symbol(:g_, i) for i = 1:N]
    calls = [:($(g[i]) = ComplexF64(gfs[$i](t))) for i = 1:N]
    entries = Expr[]
    for j = 1:N, i = 1:N  # column-major
        e = if i == j
            :(zero(ComplexF64))
        elseif i < j
            :(0.5 * $(g[i]) * conj($(g[j])))
        else
            :(-0.5 * conj($(g[j])) * $(g[i]))
        end
        push!(entries, e)
    end
    quote
        $(calls...)
        SMatrix{$N,$N,ComplexF64,$(N * N)}(($(entries...),))
    end
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
