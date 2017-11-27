# BitIndex
# --------
#
# Utils for indexing bits in a vector of 64-bit integers (internal use only).
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

struct BitIndex
    val::Int64
end

BitIndex(index, nbits) = BitIndex((index - 1) << trailing_zeros(nbits))

#         index(i)-1        index(i)        index(i)+1
# ....|................|..X.............|................|....
#                          |<-offset(i)-|
#                      |<--- 64 bits -->|

index(i::BitIndex) = (i.val >> 6) + 1
offset(i::BitIndex) = i.val & 0b111111

Base.:+(i::BitIndex, n::Int) = BitIndex(i.val + n)
Base.:-(i::BitIndex, n::Int) = BitIndex(i.val - n)
Base.:-(i1::BitIndex, i2::BitIndex) = i1.val - i2.val
Base.:(==)(i1::BitIndex, i2::BitIndex) = i1.val == i2.val
Base.isless(i1::BitIndex, i2::BitIndex) = isless(i1.val, i2.val)
Base.cmp(i1::BitIndex, i2::BitIndex) = cmp(i1.val, i2.val)

function Base.iterate(i::BitIndex, s=1)
    if s == 1
        return (index(i), 2)
    elseif s == 2
        return (offset(i), 3)
    else
        return nothing
    end
end

Base.show(io::IO, i::BitIndex) = print(io, '(', index(i), ", ", offset(i), ')')

# Create a bit mask that fills least significant `n` bits (`n` must be a
# non-negative integer).
bitmask(::Type{A}) where {A <: Alphabet} = bitmask(bitsof(A))
bitmask(n::Integer) = bitmask(UInt64, n)
bitmask(::Type{T}, n::Integer) where {T} = (one(T) << n) - one(T)
