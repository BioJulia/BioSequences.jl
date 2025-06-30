const BitUnsigned = Union{UInt8, UInt16, UInt32, UInt64, UInt128}

@inline function reversebits(x::T, ::BitsPerSymbol{2}) where T <: BitUnsigned
     mask = 0x33333333333333333333333333333333 % T
     x = ((x >> 2) & mask) | ((x & mask) << 2)
     return reversebits(x, BitsPerSymbol{4}())
end

@inline function reversebits(x::T, ::BitsPerSymbol{4}) where T <: BitUnsigned
     mask = 0x0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F % T
     x = ((x >> 4) & mask) | ((x & mask) << 4)
     return reversebits(x, BitsPerSymbol{8}())
end

@inline reversebits(x::T, ::BitsPerSymbol{8}) where T <: BitUnsigned = bswap(x)

@inline reversebits(x::UInt16, ::BitsPerSymbol{16}) = x
@inline function reversebits(x::T, ::BitsPerSymbol{16}) where T <: Union{UInt32, UInt64}
    mask = 0x0000FFFF0000FFFF0000FFFF0000FFFF % T
    x = ((x >> 16) & mask) | ((x & mask) << 16)
    reversebits(x, BitsPerSymbol{32}())
end

@inline reversebits(x::UInt32, ::BitsPerSymbol{32}) = x
@inline function reversebits(x::T, ::BitsPerSymbol{32}) where T <: Union{UInt64}
    mask = 0x00000000FFFFFFF00000000FFFFFFFF % T
    x = ((x >> 32) & mask) | ((x & mask) << 32)
    reversebits(x, BitsPerSymbol{64}())
end

@inline reversebits(x::UInt64, ::BitsPerSymbol{64}) = x

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
    return count_ones((x ⊻ (x >>> 1)) & msk)
end

@inline function gc_bitcount(x::Unsigned, ::NucleicAcidAlphabet{4})
    a =  x
    c = x >> 1
    g = x >> 2
    t = x >> 3
    gc = (c | g) & ~(a | t)
    return count_ones(gc & repeatpattern(typeof(x), 0x11))
end

@inline function ambiguous_bitcount(x::UInt64, ::T) where {T<:NucleicAcidAlphabet{4}}
    return count_nonzero_nibbles(enumerate_nibbles(x) & 0xEEEEEEEEEEEEEEEE)
end

@inline function ambiguous_bitcount(x::UInt64, ::AminoAcidAlphabet)
    return count_compared_bytes(i -> (i > 0x15) & (i < 0x1a), x)
end

@inline function count_compared_bytes(f::F, x::UInt64) where F
    z = @inline reinterpret(NTuple{8, VecElement{UInt8}}, x)
    @inline sum(map(i -> f(i.value), z))
end

pattern(::BitsPerSymbol{1})   = typemax(UInt128)
pattern(::BitsPerSymbol{2})   = 0x55555555555555555555555555555555
pattern(::BitsPerSymbol{4})   = 0x11111111111111111111111111111111
pattern(::BitsPerSymbol{8})   = 0x01010101010101010101010101010101
pattern(::BitsPerSymbol{16})  = 0x00010001000100010001000100010001
pattern(::BitsPerSymbol{32})  = 0x00000001000000010000000100000001
pattern(::BitsPerSymbol{64})  = 0x00000000000000010000000000000001
pattern(::BitsPerSymbol{128}) = 0x00000000000000000000000000000001

function count_encoding(chunk::T, encoding::T, b::BitsPerSymbol{B}) where {T <: Unsigned, B}
    pat = pattern(b) % typeof(encoding)
    u = chunk ⊻ (encoding * pat)
    for i in 1:trailing_zeros(B)
        shift = (1 << (i - 1)) & (8*sizeof(T) - 1)
        u |= (u >> shift)
    end
    u = ~u & pat
    return count_ones(u)
end

@inline function certain_bitcount(x::UInt64, ::T) where {T<:NucleicAcidAlphabet{4}}
    x = enumerate_nibbles(x)
    x = x ⊻ 0x1111111111111111
    return count_0000_nibbles(x)
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
