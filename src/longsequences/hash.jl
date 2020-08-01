###
### Hash
###
###
### MurmurHash3 function of BioSequence.
###
### The hash function defined here cares about the starting position of a
### character sequence in the underlying data. That means, even if the starting
### positions of two sequences (`s1` and `s2`) are different in their `data`
### field, their hash values are identical if `s1 == s1` is true.  The
### implementation is based on the 128bit MurmurHash3 function, which was written
### by Austin Appleby, and the source code is distributed under the public domain:
### https://github.com/aappleby/smhasher/blob/master/src/MurmurHash3.cpp
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

const c1 = 0x87c37b91114253d5
const c2 = 0x4cf5ad432745937f

@inline function rotl64(x::UInt64, r)
    return (x << (r & 63)) | (x >>> (-r & 63))
end

@inline function fmix64(k::UInt64)
    k = k ⊻ k >> 33
    k *= 0xff51afd7ed558ccd
    k = k ⊻ k >> 33
    k *= 0xc4ceb9fe1a85ec53
    k = k ⊻ k >> 33
    return k
end

@inline function murmur1(h1, k1)
    k1 *= c1
    k1 = rotl64(k1, 31)
    k1 *= c2
    h1 = h1 ⊻ k1
    return (h1, k1)
end

@inline function murmur2(h1, h2, k2)
    k2 *= c2
    k2 = rotl64(k2, 33)
    k2 *= c1
    h2 = h1 ⊻ k2
    return (h2, k2)
end

@inline function murmur(h1, h2, k1, k2)
    h1, k1 = murmur1(h1, k1)
    h1 = rotl64(h1, 27)
    h1 += h2
    h1 = h1 * 5 + 0x52dce729

    h2, k2 = murmur2(h1, h2, k2)
    h2 = rotl64(h2, 31)
    h2 += h1
    h2 = h2 * 5 + 0x38495ab5

    return (h1, h2, k1, k2)
end

function finalize(h1, h2, len)
    h1 = h1 ⊻ len
    h2 = h2 ⊻ len
    h1 += h2
    h2 += h1
    h1 = fmix64(h1)
    h2 = fmix64(h2)
    h1 += h2
    h2 += h1

    return h1
end

function tail(data, next, stop, h1, h2)
    k1 = k2 = zero(UInt64)
    
    # Use this to mask any noncoding bits in last chunk
    mask = bitmask(offset(stop - bits_per_symbol(stop)) + bits_per_symbol(stop))
    @inbounds if next < stop
        k1 = data[index(next)]
        next += 64
    end
    @inbounds if next < stop
        k2 = data[index(next)] & mask
    else
        k1 &= mask
    end
    
    h1, k1 = murmur1(h1, k1)
    h2, k2 = murmur2(h1, h2, k2)
    return (h1, h2)
end

# ref: MurmurHash3_x64_128
function Base.hash(seq::LongSequence, seed::UInt64)
    # Mix sequence length so that dna"A" and dna"AA"
    # return the different hash values.
    h1::UInt64 = h2::UInt64 = hash(length(seq), seed)
    next = bitindex(seq, 1)
    last = bitindex(seq, lastindex(seq) + 1)
    data = encoded_data(seq)

    @inbounds while last - next ≥ 128
        j = index(next)
        k1 = data[j]
        k2 = data[j+1]
        h1, h2, k1, k2 = murmur(h1, h2, k1, k2)
        next += 128
    end

    h1, h2 = tail(data, next, last, h1, h2)
    h1 = finalize(h1, h2, length(seq))

    return h1
end
