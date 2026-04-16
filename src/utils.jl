"""
    substitute_operators(op, dict::Dict; replace_adjoint=true)

Like `substitute(op, dict::Dict)` but with special handling for `QMul` and `QAdd`.
This is needed if an operator is substitute by a `QMul` or `QAdd`, e.g.
    ``a_1 \\rightarrow g_2 a_2 + g_3 a_3``

If `replace_adjoint=true`, the dictionary is extended with adjoint substitutions
for all key/value pairs, i.e. `adjoint(key) => adjoint(value)`.
"""
function substitute_operators(op, dict::Dict; replace_adjoint = true)
    dict_ = replace_adjoint ? _extend_with_adjoint(dict) : dict
    return substitute(op, dict_)
end
function substitute_operators(op::SQA.QAdd, dict::Dict; replace_adjoint = true)
    dict_ = replace_adjoint ? _extend_with_adjoint(dict) : dict
    return SQA.QAdd([
        substitute_operators(arg, dict_; replace_adjoint = false) for arg in op.arguments
    ])
end
function substitute_operators(op::SQA.QMul, dict::Dict; replace_adjoint = true)
    dict_ = replace_adjoint ? _extend_with_adjoint(dict) : dict
    return substitute(op.arg_c, dict_) * prod([substitute(arg, dict_) for arg ∈ op.args_nc])
end

function _extend_with_adjoint(dict::Dict)
    pairs_ = collect(dict)
    adj_pairs = [Base.adjoint(k) => Base.adjoint(v) for (k, v) in pairs_]
    return Dict(vcat(pairs_, adj_pairs))
end
