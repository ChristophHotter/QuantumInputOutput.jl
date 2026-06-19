"""
    substitute_operators(op, dict::Dict; replace_adjoint=true)

Like `substitute(op, dict::Dict)` but with special handling for operator products and sums.
This is needed if an operator is substituted by a product or sum of operators, e.g.
    ``a_1 \\rightarrow g_2 a_2 + g_3 a_3``

Plain `substitute` recurses into the replacement, so substituting an operator by an
expression that still contains that operator loops forever. This routine instead replaces
each leaf operator with a single dictionary lookup and rebuilds the product/sum.

If `replace_adjoint=true`, the dictionary is extended with adjoint substitutions
for all key/value pairs, i.e. `adjoint(key) => adjoint(value)`.
"""
function substitute_operators(op, dict::Dict; replace_adjoint = true)
    dict_ = replace_adjoint ? _extend_with_adjoint(dict) : dict
    return substitute(op, dict_)
end

function substitute_operators(op::SQA.QSym, dict::Dict; replace_adjoint = true)
    dict_ = replace_adjoint ? _extend_with_adjoint(dict) : dict
    return get(dict_, op, op)
end

function substitute_operators(op::SQA.QAdd, dict::Dict; replace_adjoint = true)
    dict_ = replace_adjoint ? _extend_with_adjoint(dict) : dict
    iszero(op) && return op
    return sum(c * _replace_ops(term.ops, dict_) for (term, c) in op)
end

# Single-level operator replacement: each leaf is replaced by a direct dictionary
# lookup (never re-substituted), avoiding the infinite recursion of `substitute`.
_replace_ops(ops, dict_) = isempty(ops) ? 1 : prod(get(dict_, o, o) for o in ops)

function _extend_with_adjoint(dict::Dict)
    pairs_ = collect(dict)
    adj_pairs = [Base.adjoint(k) => Base.adjoint(v) for (k, v) in pairs_]
    return Dict(vcat(pairs_, adj_pairs))
end
