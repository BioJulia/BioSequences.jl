
# Counters
# --------

count_gc(x::ShortSequence) = gc_bitcount(encoded_data(x), BitsPerSymbol(x))

function count_a(x::ShortSequence)
    # Have to take into account 00 bitpairs that are part of the unused bit of
    # the uint. Hence the subtraction.
    return count_a(encoded_data(x)) - n_unused(x)
end

function count_c(x::ShortSequence)
    return count_c(encoded_data(x))
end

function count_g(x::ShortSequence)
    return count_g(encoded_data(x))
end

function count_t(x::ShortSequence)
    return count_t(encoded_data(x))
end

"""
    mismatches(a::ShortSequence, b::ShortSequence)

Return the number of mismatches between `a` and `b`.
"""
function mismatches(a::ShortSequence, b::ShortSequence)
    return count_nonzero_bitpairs(encoded_data(a) ⊻ encoded_data(b))
end