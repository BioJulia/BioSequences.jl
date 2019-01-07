# Transformations
# ---------------

"""
    complement(x::T) where {T <: Skipmer}

Return the complement of a short sequence type `x`.
"""
function BioSymbols.complement(x::T) where {T <: Skipmer}
    return T(~encoded_data(x))
end


"""
    reverse(x::T) where {T <: Skipmer}

Return the reverse of short sequence type variable `x`.
"""
function Base.reverse(x::T) where {T <: Skipmer}
    bits = encoded_data(x)
    rbits = reversebits(bits, BitsPerSymbol{2}())
    return T(rbits >> (sizeof(bits) * 8 - 2 * length(x)))
end


"""
    reverse_complement(x::Skipmer)

Return the reverse complement of `x`.
"""
reverse_complement(x::Skipmer) = complement(reverse(x))


"""
    canonical(kmer::Skipmer)

Return the canonical sequence of `x`.

A canonical sequence is the numerical lesser of a k-mer and its reverse complement.
This is useful in hashing/counting sequences in data that is not strand specific,
and thus observing the short sequence is equivalent to observing its reverse complement.
"""
function canonical(x::Skipmer)
    y = reverse_complement(x)
    return x < y ? x : y
end


function swap(x::T, i, j) where {T <: Skipmer}
    i = 2 * length(x) - 2i
    j = 2 * length(x) - 2j
    b = encoded_data(x)
    x = ((b >> i) ⊻ (b >> j)) & encoded_data_eltype(x)(0x03)
    return T(b ⊻ ((x << i) | (x << j)))
end


function Random.shuffle(x::T) where {T <: Skipmer}
    # Fisher-Yates shuffle
    j = lastindex(x)
    for i in firstindex(x):(j - 1)
        j′ = rand(i:j)
        x = swap(x, i, j′)
    end
    return x
end