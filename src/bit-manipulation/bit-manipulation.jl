
include("bitindex.jl")
include("bitpar-compiler.jl")

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
