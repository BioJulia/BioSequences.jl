###
### Mer specific specializations of src/biosequence/transformations.jl
###

"""
    complement(x::T) where {T <: Skipmer}

Return the complement of a short sequence type `x`.
"""
function BioSymbols.complement(x::T) where {T<:AbstractMer}
    return reinterpret(T, complement_bitpar(encoded_data(x), Alphabet(x)))
end

"""
    reverse(x::T) where {T <: Skipmer}

Return the reverse of short sequence type variable `x`.
"""
function Base.reverse(x::T) where {T<:AbstractMer}
    bits = encoded_data(x)
    rbits = reversebits(bits, BitsPerSymbol{2}())
    return reinterpret(T, rbits >> (sizeof(bits) * 8 - 2 * length(x)))
end


"""
    reverse_complement(x::Skipmer)

Return the reverse complement of `x`.
"""
reverse_complement(x::AbstractMer) = complement(reverse(x))


"""
    canonical(kmer::Skipmer)

Return the canonical sequence of `x`.

A canonical sequence is the numerical lesser of a k-mer and its reverse complement.
This is useful in hashing/counting sequences in data that is not strand specific,
and thus observing the short sequence is equivalent to observing its reverse complement.
"""
@inline canonical(x::AbstractMer) = min(x, reverse_complement(x))


function swap(x::T, i, j) where {T<:AbstractMer}
    i = 2 * length(x) - 2i
    j = 2 * length(x) - 2j
    b = encoded_data(x)
    x = ((b >> i) ⊻ (b >> j)) & encoded_data_type(x)(0x03)
    return reinterpret(T, b ⊻ ((x << i) | (x << j)))
end


function Random.shuffle(x::T) where {T<:AbstractMer}
    # Fisher-Yates shuffle for mers.
    j = lastindex(x)
    for i in firstindex(x):(j - 1)
        j′ = rand(i:j)
        x = swap(x, i, j′)
    end
    return x
end
