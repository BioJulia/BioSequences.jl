###
### Mer: Predicates & comparisons
###

Base.cmp(seq1::T, seq2::T) where {T<:AbstractMer} = cmp(encoded_data(seq1), encoded_data(seq2))
Base.:(==)(x::T, y::T) where {T<:AbstractMer} = encoded_data(x) == encoded_data(y)
Base.isless(x::T, y::T) where {T<:AbstractMer} = isless(encoded_data(x), encoded_data(y))
