###
### Mer specific specializations of src/biosequence/predicates.jl
###

Base.cmp(seq1::T, seq2::T) where {T<:AbstractMer} = cmp(encoded_data(seq1), encoded_data(seq2))
Base.:(==)(x::T, y::T) where {T<:AbstractMer} = encoded_data(x) == encoded_data(y)
Base.isless(x::T, y::T) where {T<:AbstractMer} = isless(encoded_data(x), encoded_data(y))

function Base.hash(x::Mer{<:NucleicAcidAlphabet{2},K}, h::UInt) where {K}
    return Base.hash_uint64(encoded_data(x) ⊻ K ⊻ h)
end
