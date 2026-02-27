
# use tuple (x,y) instead of vector [x,y] if the type is important (L?)

"""
    SLH

SLH triple with scattering matrix ``S``, Lindblad term ``L`` and Hamiltonian ``H``. 
``S`` and ``L`` can also be vectors of scattering matrices and Linblad terms

See also [`▷`](@ref) and [`⊞`](@ref)
"""
struct SLH{T,LT,H,S<:AbstractMatrix{T},L<:AbstractVector{LT}} # TODO: L[i] and H QTerm or Number
    scattering::S
    lindblad::L
    hamiltonian::H

    function SLH(scattering::S, lindblad::L, hamiltonian::H) where
        {T,LT,H,S<:AbstractMatrix{T},L<:AbstractVector{LT}}
        
        @assert size(scattering, 1) == length(lindblad)
        new{T,LT,H,S,L}(scattering, lindblad, hamiltonian)
    end
end

# SLH(scattering::S, lindblad::L, hamiltonian::H) where {S,H,L} = SLH{S,H,L}(scattering,lindblad,hamiltonian)
function SLH(S,L::AbstractVector,H) 
    S_ = Matrix{typeof(S)}(I,length(L),length(L))*S 
    return SLH(S_,L,H)
end
function SLH(S,L,H) 
    S_ = ones(1,1)*S
    L_ = [L]
    return SLH(S_,L_,H)
end

function Base.isequal(slh1::SLH, slh2::SLH)
    isequal(get_scattering(slh1), get_scattering(slh2)) && isequal(get_lindblad(slh1), get_lindblad(slh2)) && isequal(get_hamiltonian(slh1), get_hamiltonian(slh2))
end

get_scattering(slh::SLH) = slh.scattering
get_lindblad(slh::SLH) = slh.lindblad
get_hamiltonian(slh::SLH) = slh.hamiltonian

"""
    ▷(G::SLH...)

Creates a new SLH triple by cascading the SLH triples from first to last according to 
the rule ``SLH_1 \\triangleright SLH_2 = ( S_2 S1, L_2 + S_2 L_1, H_1 + H_2 - \\frac{i}{2} L_2^\\dagger S_2 L_1 - L_1^\\dagger S_2^\\dagger L_2 ) ``
Unicode `\\triangleright<tab>` alias of [`cascade`](@ref)
"""
function ▷(G1::SLH,G2::SLH) #\triangleright
    # function cascade(G1::SLH,G2::SLH) # G1 ▷ G2 = G2 ◁ G1
    S1 = get_scattering(G1); L1 = get_lindblad(G1); H1 = get_hamiltonian(G1)
    S2 = get_scattering(G2); L2 = get_lindblad(G2); H2 = get_hamiltonian(G2)
    
    lL1 = length(L1)
    lL2 = length(L2)
    @assert lL1==lL2

    S_t = simplify.(S2*S1)
    L_t = simplify.(L2 + S2*L1)
    # H_t = (H1 + H2 - 1im/2*(L2'S2*L1 - L1'S2'L2)) # the adjoint also creates the transpose! 
    # H_t = (H1 + H2 - 1im/2*(QC._adjoint(L2)*S2*L1 - QC._adjoint(L1)*QC._adjoint(S2)*L2))
    L1_ct = QC._adjoint.(L1)
    L2_ct = QC._adjoint.(L2)
    S2_ct = QC._adjoint.(S2)

    H_t = simplify(H1 + H2 - 1im/2*( sum(L2_ct[it]*(S2*L1)[it] for it=1:lL1 ) - 
            sum(L1_ct[it]*(S2_ct*L2)[it] for it=1:lL1 ) ))
    return SLH(S_t, L_t, H_t)
end
▷(a::SLH, b::SLH, c::SLH...) = ▷(a▷b,c...)
# G's in the order as they are cascaded (from left to right)

# TODO: vectors S, L
"""
    cascade(G::SLH...)

Creates a new SLH triple by cascading the SLH triples from first to last according to 
the rule ``SLH_1 \\triangleright SLH_2 = ( S_2 S1, L_2 + S_2 L_1, H_1 + H_2 - \\frac{i}{2} L_2^\\dagger S_2 L_1 - L_1^\\dagger S_2^\\dagger L_2 ) ``
See also [`▷`](@ref). 
"""
cascade(args...) = ▷(args...)
# ◁(G1::SLH,G2::SLH) = ▷(G2,G1) # unkwnown unicode character

# TODO: equation for concateneate
"""
    ⊞(G::SLH...)

Creates a new SLH triple by concatenating the SLH triples according to 
the rule ``SLH_1 \\boxplus SLH_2 = TODO``
Unicode `\\boxplus<tab>` alias of [`concatenate`](@ref)
"""
function ⊞(G1::SLH,G2::SLH) #\boxplus
    S1 = get_scattering(G1); L1 = get_lindblad(G1); H1 = get_hamiltonian(G1)
    S2 = get_scattering(G2); L2 = get_lindblad(G2); H2 = get_hamiltonian(G2)

    lS1 = size(S1)[1]; lS2 = size(S2)[1]

    T_S = promote_type(eltype(S1), eltype(S2))
    S_t = zeros(T_S, lS1+lS2, lS1+lS2)
    S_t[1:lS1,1:lS1] .= S1
    S_t[1+lS1:lS1+lS2,1+lS1:lS1+lS2] .= S2

    L_t = [L1; L2]
    H_t = H1 + H2

    return SLH(S_t, L_t, H_t)
end
⊞(a::SLH, b::SLH, c::SLH...) = ⊞(a⊞b,c...)

"""
    concatenate(G::SLH...)

Creates a new SLH triple by concatenating the SLH triples according to 
the rule ``SLH_1 \\boxplus SLH_2 = TODO``
See also [`⊞`](@ref).
"""
concatenate(args...) = ⊞(args...)
concatenation(args...) = concatenate(args...) # alias

Base.length(h::SecondQuantizedAlgebra.ConcreteHilbertSpace) = 1
Base.length(h::ProductSpace) = length(h.spaces)


### direct QuantumOptics implementation of SLH to QuantumOptics.jl operators
### translate SLH to (time-dependent) Hamiltonian and Lindblad operators

# struct SLHqo{S<:AbstractMatrix,L<:AbstractVector,H}
#     scattering::S
#     lindblad::L
#     hamiltonian::H
# end
# function SLHqo(scattering::S, lindblad::L, hamiltonian::H) where {S<:AbstractMatrix,L<:AbstractVector,H}
#     @assert size(scattering, 1) == length(lindblad)
#     return SLHqo{S,L,H}(scattering, lindblad, hamiltonian)
# end

# dispatch SLH(S,L,H::Op or function) to SLHqo
