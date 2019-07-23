
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

@inline function complementbits(x::UInt64, ::T) where {T <: NucleicAcidAlphabet{2}}
    return ~x
end

@inline function complementbits(x::UInt64, ::T) where {T <: NucleicAcidAlphabet{4}}
    return (
        ((x & 0x1111111111111111) << 3) | ((x & 0x8888888888888888) >> 3) |
        ((x & 0x2222222222222222) << 1) | ((x & 0x4444444444444444) >> 1)
    )
end

@inline function gc_bitcount(x::UInt64, ::BitsPerSymbol{2})
    msk = 0x5555555555555555
    c = x & msk
    g = (x >> 1) & msk
    return count_ones(c ⊻ g)
end

@inline function gc_bitcount(x::UInt64, ::BitsPerSymbol{4})
    a =  x & 0x1111111111111111
    c = (x & 0x2222222222222222) >> 1
    g = (x & 0x4444444444444444) >> 2
    t = (x & 0x8888888888888888) >> 3
    return count_ones((c | g) & ~(a | t))
end

@inline count_a(x::UInt64) = count_00_bitpairs(x)
@inline count_c(x::UInt64) = count_01_bitpairs(x)
@inline count_g(x::UInt64) = count_10_bitpairs(x)
@inline count_t(x::UInt64) = count_11_bitpairs(x)