
include("bitindex.jl")
include("bitpar-compiler.jl")

@inline function reversebits(x::UInt64, ::BitsPerSymbol{2})
     x =  ((x >> 2) & 0x3333333333333333) | ((x & 0x3333333333333333) << 2)
     x = reversebits(x, BitsPerSymbol{4}())
     return x
end

@inline function reversebits(x::UInt64, ::BitsPerSymbol{4})
     x = ((x >> 4 ) & 0x0F0F0F0F0F0F0F0F) | ((x & 0x0F0F0F0F0F0F0F0F) << 4 )
     x = ((x >> 8 ) & 0x00FF00FF00FF00FF) | ((x & 0x00FF00FF00FF00FF) << 8 )
     x = ((x >> 16) & 0x0000FFFF0000FFFF) | ((x & 0x0000FFFF0000FFFF) << 16)
     x = ((x >> 32) & 0x00000000FFFFFFFF) | ((x & 0x00000000FFFFFFFF) << 32)
     return x
end

@inline function reversebits(x::UInt128, ::BitsPerSymbol{2})
     x =  ((x >> 2) & 0x33333333333333333333333333333333) | ((x & 0x33333333333333333333333333333333) << 2)
     x = reversebits(x, BitsPerSymbol{4}())
     return x
end

@inline function reversebits(x::UInt128, ::BitsPerSymbol{4})
     x = ((x >> 4 ) & 0x0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F) | ((x & 0x0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F) << 4 )
     x = ((x >> 8 ) & 0x00FF00FF00FF00FF00FF00FF00FF00FF) | ((x & 0x00FF00FF00FF00FF00FF00FF00FF00FF) << 8 )
     x = ((x >> 16) & 0x0000FFFF0000FFFF0000FFFF0000FFFF) | ((x & 0x0000FFFF0000FFFF0000FFFF0000FFFF) << 16)
     x = ((x >> 32) & 0x00000000FFFFFFFF00000000FFFFFFFF) | ((x & 0x00000000FFFFFFFF00000000FFFFFFFF) << 32)
     x = ((x >> 64) & 0x0000000000000000FFFFFFFFFFFFFFFF) | ((x & 0x0000000000000000FFFFFFFFFFFFFFFF) << 64)
     return x
end

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

@inline count_a(x::Unsigned) = count_00_bitpairs(x)
@inline count_c(x::Unsigned) = count_01_bitpairs(x)
@inline count_g(x::Unsigned) = count_10_bitpairs(x)
@inline count_t(x::Unsigned) = count_11_bitpairs(x)