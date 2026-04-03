###########################
### interaction picture ###
###########################

using StaticArrays: SMatrix, MMatrix

# [I4 fix]: dispatch-based instead of applicable()
_as_time_function(x::Number) = _ -> x
_as_time_function(x) = x  # anything callable passes through

"""
    coupling_matrix(gs::Tuple)

Build the antisymmetric coupling coefficient matrix `A(t)` from a tuple of
coupling functions/constants `gs`. Returns a closure `t -> SMatrix{N,N,ComplexF64}`.

Uses `MMatrix` internally for zero-allocation construction, returns `SMatrix`.
"""
function coupling_matrix(gs::NTuple{N}) where {N}  # [I6 fix]: NTuple{N} not NTuple{N,Any}
    gfs = map(_as_time_function, gs)
    function A(t)
        M = zero(MMatrix{N,N,ComplexF64})
        @inbounds for i in 1:N, j in (i+1):N
            val = 0.5 * gfs[i](t) * conj(gfs[j](t))
            M[i, j] = val
            M[j, i] = -conj(val)
        end
        return SMatrix(M)
    end
    return A
end

coupling_matrix(g1, g2, gs...) = coupling_matrix((g1, g2, gs...))

"""
    solve_mode_evolution(A::Function, T; alg=Tsit5(), kwargs...)

Solve the interaction-picture coefficient-matrix ODE `dM/dt = A(t) M(t)` with `M(0) = I`.
Uses in-place ODE formulation (`mul!`) for zero allocations per step.

Returns the ODE solution directly (callable as `sol(t)`).
"""
function solve_mode_evolution(A::Function, T; alg = Tsit5(), kwargs...)
    T0 = T[1]
    Tend = T[end]
    n = size(A(T0), 1)
    M0 = Matrix{ComplexF64}(I, n, n)
    # [I3 fix]: no _A_buf, mul! works directly with SMatrix
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
Returns `t -> SMatrix{2,2,Float64}`.
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
