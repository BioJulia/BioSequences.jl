###
### Mer: Counters
###

count_gc(x::AbstractMer) = gc_bitcount(encoded_data(x), BitsPerSymbol{2}())

function count_a(x::AbstractMer)
    # Have to take into account 00 bitpairs that are part of the unused bit of
    # the uint. Hence the subtraction.
    return count_a(encoded_data(x)) - n_unused(x)
end

function count_c(x::AbstractMer)
    return count_c(encoded_data(x))
end

function count_g(x::AbstractMer)
    return count_g(encoded_data(x))
end

function count_t(x::AbstractMer)
    return count_t(encoded_data(x))
end

"""
    mismatches(a::AbstractMer, b::AbstractMer)

Return the number of mismatches between `a` and `b`.
"""
function mismatches(a::T, b::T) where {T<:AbstractMer}
    return count_nonzero_bitpairs(encoded_data(a) ⊻ encoded_data(b))
end
