###########################
### interaction picture ###
###########################

### get matrix M(t)
function get_Mt(A::Function, T)
    T0 = T[1]
    Tend = T[end]
    l_A = size(A(T0))[1]
    M0 = diagm(ones(ComplexF64, l_A))
    f_M(u, p, t) = A(t) * u
    prob_M = ODEProblem(f_M, M0, (T0, Tend))
    return solve(prob_M; saveat = T) #, abstol, reltol) # TODO: kwargs
end

# ### get matrix A(t) (for 2x2, 3x3 and 4x4 matrices)
function get_At_2(gu, gv) # assumes coupling u-v
    A(t_) = 0.5 * [0 gu(t_) * gv(t_)';
        -gu(t_)' * gv(t_) 0]
    return A
end
get_At_2(g_ls) = get_At_2(g_ls...)

function get_At_3(g_ls::Vector) # assumes coupling u-s-v
    @assert length(g_ls) == 3
    g1 = g_ls[1]
    g2 = g_ls[2]
    g3 = g_ls[3]
    A(t_) = 0.5 * [0 g1(t_) * g2(t_)' g1(t_) * g3(t_)';
        -g1(t_)' * g2(t_) 0 g2(t_) * g3(t_)';
        -g1(t_)' * g3(t_) -g2(t_)' * g3(t_) 0]
    return A
end

# gu_ls: list of gu functions for the corresponding operators au
# gv_ls: list of gv functions for the corresponding operators av
function get_At_4(gu_ls::Vector, gv_ls::Vector)
    @assert length(gu_ls) + length(gv_ls) == 4
    gu1 = gu_ls[1]
    gu2 = gu_ls[2]
    gv1 = gv_ls[1]
    gv2 = gv_ls[2]
    A(t_) = 0.5 * [0 -gu1(t_) * gu2(t_)' gu1(t_) * gv1(t_)' gu1(t_) * gv2(t_)';
        gu1(t_)' * gu2(t_) 0 gu2(t_) * gv1(t_)' gu2(t_) * gv2(t_)';
        -gu1(t_)' * gv1(t_) -gu2(t_)' * gv1(t_) 0 gv1(t_) * gv2(t_)';
        -gu1(t_)' * gv2(t_) -gu2(t_)' * gv2(t_) -gv1(t_)' * gv2(t_) 0]
    return A
end
