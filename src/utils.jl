"""
    substitute_operators(op, dict::Dict; replace_adjoint=true)

Like `substitute(op, dict::Dict)`, but extends the substitution dictionary with
adjoint mappings when requested. This is useful when an operator replacement
contains sums or products, e.g. ``a_1 \\rightarrow g_2 a_2 + g_3 a_3``.

If `replace_adjoint=true`, the dictionary is extended with adjoint substitutions
for all key/value pairs, i.e. `adjoint(key) => adjoint(value)`.
"""
function substitute_operators(op, dict::Dict; replace_adjoint = true)
    dict_ = replace_adjoint ? _extend_with_adjoint(dict) : dict
    return substitute(op, dict_)
end

function substitute_operators(op::SQA.QAdd, dict::Dict; replace_adjoint = true)
    dict_ = replace_adjoint ? _extend_with_adjoint(dict) : dict
    result = zero(SQA.QAdd)
    for (term, c) in op.arguments
        term_expr = c
        for arg in term.ops
            term_expr *= get(dict_, arg, arg)
        end
        result += term_expr
    end
    return result
end

function _extend_with_adjoint(dict::Dict)
    pairs_ = collect(dict)
    adj_pairs = [Base.adjoint(k) => Base.adjoint(v) for (k, v) in pairs_]
    return Dict(vcat(pairs_, adj_pairs))
end
