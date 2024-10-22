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

@inline function gc_bitcount(x::Unsigned, ::NucleicAcidAlphabet{2})
    msk = repeatpattern(typeof(x), 0x55)
    c = x & msk
    g = (x >> 1) & msk
    return count_ones(c ⊻ g)
end

@inline function gc_bitcount(x::Unsigned, ::NucleicAcidAlphabet{4})
    a =  x & repeatpattern(typeof(x), 0x11)
    c = (x & repeatpattern(typeof(x), 0x22)) >> 1
    g = (x & repeatpattern(typeof(x), 0x44)) >> 2
    t = (x & repeatpattern(typeof(x), 0x88)) >> 3
    return count_ones((c | g) & ~(a | t))
end

@inline a_bitcount(x::Unsigned, ::NucleicAcidAlphabet{2}) = count_00_bitpairs(x)
@inline c_bitcount(x::Unsigned, ::NucleicAcidAlphabet{2}) = count_01_bitpairs(x)
@inline g_bitcount(x::Unsigned, ::NucleicAcidAlphabet{2}) = count_10_bitpairs(x)
@inline t_bitcount(x::Unsigned, ::NucleicAcidAlphabet{2}) = count_11_bitpairs(x)

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

@inline function ambiguous_bitcount(x::UInt64, ::AminoAcidAlphabet)
    return count_compared_bytes(i -> (i > 0x15) & (i < 0x1a), x)
end

@inline function ambiguous_bitcount(a::UInt64, b::UInt64, ::T) where {T<:NucleicAcidAlphabet{4}}
    return count_nonzero_nibbles((enumerate_nibbles(a) | enumerate_nibbles(b)) & 0xEEEEEEEEEEEEEEEE)
end

@inline function gap_bitcount(x::UInt64, ::T) where {T<:NucleicAcidAlphabet{4}}
    return count_0000_nibbles(x)
end

@inline function gap_bitcount(a::UInt64, b::UInt64, ::T) where {T<:NucleicAcidAlphabet{4}}
    # Count the gaps in a, count the gaps in b, subtract the number of shared gaps.
    return count_0000_nibbles(a) + count_0000_nibbles(b) - count_0000_nibbles(a | b)
end

@inline function count_compared_bytes(f::F, x::UInt64) where F
    z = @inline reinterpret(NTuple{8, VecElement{UInt8}}, x)
    @inline sum(map(i -> f(i.value), z))
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

@inline function four_to_two_bits(x::UInt64)::UInt32
    m1 = 0x1111111111111111
    m2 = m1 << 1
    m3 = m1 | m2
    y  = (x >> 3) & m1
    y |= (x >> 2) & m2
    y |= (x >> 1) & m3
    pack(y)
end

@inline function pack(x::UInt64)::UInt32
    m1 = 0x0f0f0f0f0f0f0f0f
    m2 = 0x00ff00ff00ff00ff
    m3 = 0x0000ffff0000ffff
    m4 = 0x00000000ffffffff
    x = (x & m1) | (x & ~m1) >> 2
    x = (x & m2) | (x & ~m2) >> 4
    x = (x & m3) | (x & ~m3) >> 8
    x = (x & m4) | (x & ~m4) >> 16
    x % UInt32
end

@inline function two_to_four_bits(x::UInt32)::UInt64
    m = 0x1111111111111111
    y = expand(x)
    z = (y & m) << 1 | (m & ~y)
    m2 = y & (m << 1)
    m2 = m2 | m2 >> 1
    (z & m2) << 2 | (z & ~m2)
end

@inline function expand(x::UInt32)::UInt64
    m1 = 0x000000000000ffff
    m2 = 0x000000ff000000ff
    m3 = 0x000f000f000f000f
    m4 = 0x0303030303030303
    y = x % UInt64
    y = (y & m1) | (y & ~m1) << 16
    y = (y & m2) | (y & ~m2) << 8
    y = (y & m3) | (y & ~m3) << 4
        (y & m4) | (y & ~m4) << 2
end
