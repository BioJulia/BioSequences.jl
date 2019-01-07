# Predicates & comparisons
# ------------------------

Base.cmp(seq1::T, seq2::T) where {T <: Skipmer} = cmp(encoded_data(seq1), encoded_data(seq2))
Base.:(==)(x::T, y::T) where {T <: Skipmer} = encoded_data(x) == encoded_data(y)
Base.isless(x::T, y::T) where {T <: Skipmer} = isless(encoded_data(x), encoded_data(y))