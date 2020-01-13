
include("bitindex.jl")
include("bitpar-compiler.jl")

@inline function reversebits(x::T, ::BitsPerSymbol{2}) where T <: Base.BitUnsigned
     mask = 0x33333333333333333333333333333333 % T
     x = ((x >> 2) & mask) | ((x & mask) << 2)
     return reversebits(x, BitsPerSymbol{4}())
end

@inline function reversebits(x::T, ::BitsPerSymbol{4}) where T <: Base.BitUnsigned
     mask = 0x0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F % T
     x = ((x >> 4) & mask) | ((x & mask) << 4)
     return bswap(x)
end

reversebits(x::T, ::BitsPerSymbol{8}) where T <: Base.BitUnsigned = bswap(x)

@inline function complement_bitpar(x::Unsigned, ::T) where {T<:NucleicAcidAlphabet{2}}
    return ~x
end

@inline function complement_bitpar(x::Unsigned, ::T) where {T<:NucleicAcidAlphabet{4}}
    return (
        ((x & repeatpattern(typeof(x), 0x11)) << 3) | ((x & repeatpattern(typeof(x), 0x88)) >> 3) |
        ((x & repeatpattern(typeof(x), 0x22)) << 1) | ((x & repeatpattern(typeof(x), 0x44)) >> 1)
    )
end

@inline function gc_bitcount(x::Unsigned, ::BitsPerSymbol{2})
    msk = repeatpattern(typeof(x), 0x55)
    c = x & msk
    g = (x >> 1) & msk
    return count_ones(c ⊻ g)
end

@inline function gc_bitcount(x::Unsigned, ::BitsPerSymbol{4})
    a =  x & repeatpattern(typeof(x), 0x11)
    c = (x & repeatpattern(typeof(x), 0x22)) >> 1
    g = (x & repeatpattern(typeof(x), 0x44)) >> 2
    t = (x & repeatpattern(typeof(x), 0x88)) >> 3
    return count_ones((c | g) & ~(a | t))
end

@inline function mismatch_bitcount(a::UInt64, b::UInt64, ::T) where {T<:NucleicAcidAlphabet{4}}
    return count_nonzero_nibbles(a ⊻ b)
end

@inline function mismatch_bitcount(a::UInt64, b::UInt64, ::T) where {T<:NucleicAcidAlphabet{2}}
    return count_nonzero_bitpairs(a ⊻ b)
end

@inline function match_bitcount(a::UInt64, b::UInt64, ::T) where {T<:NucleicAcidAlphabet{4}}
    return count_0000_nibbles(a ⊻ b)
end

@inline function match_bitcount(a::UInt64, b::UInt64, ::T) where {T<:NucleicAcidAlphabet{2}}
    return count_00_bitpairs(a ⊻ b)
end

@inline function ambiguous_bitcount(x::UInt64, ::T) where {T<:NucleicAcidAlphabet{4}}
    return count_nonzero_nibbles(enumerate_nibbles(x) & 0xEEEEEEEEEEEEEEEE)
end

@inline function ambiguous_bitcount(a::UInt64, b::UInt64, ::T) where {T<:NucleicAcidAlphabet{4}}
    return count_nonzero_nibbles((enumerate_nibbles(a) | enumerate_nibbles(b)) & 0xEEEEEEEEEEEEEEEE)
end

@inline function gap_bitcount(x::UInt64, ::T) where {T<:NucleicAcidAlphabet{4}}
    return count_zero_nibbles(x)
end

@inline function gap_bitcount(a::UInt64, b::UInt64, ::T) where {T<:NucleicAcidAlphabet{4}}
    # Count the gaps in a, count the gaps in b, subtract the number of shared gaps.
    return count_0000_nibbles(a) + count_0000_nibbles(b) - count_0000_nibbles(a | b)
end

@inline function certain_bitcount(x::UInt64, ::T) where {T<:NucleicAcidAlphabet{4}}
    x = enumerate_nibbles(x)
    x = x ⊻ 0x1111111111111111
    return count_0000_nibbles(x)
end

@inline function certain_bitcount(a::UInt64, b::UInt64, ::T) where {T<:NucleicAcidAlphabet{4}}
    x = enumerate_nibbles(a) ⊻ 0x1111111111111111
    y = enumerate_nibbles(b) ⊻ 0x1111111111111111
    return count_0000_nibbles(x | y)
end

@inline count_a(x::Unsigned) = count_00_bitpairs(x)
@inline count_c(x::Unsigned) = count_01_bitpairs(x)
@inline count_g(x::Unsigned) = count_10_bitpairs(x)
@inline count_t(x::Unsigned) = count_11_bitpairs(x)
