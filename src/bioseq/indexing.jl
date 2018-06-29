# Indexing
# ========
#
# Indexing methods for biological sequences.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

# Getting and setting elements in a biological sequence.

function Base.checkbounds(seq::BioSequence, range::UnitRange)
    if 1 ≤ range.start && range.stop ≤ lastindex(seq)
        return true
    end
    throw(BoundsError(seq, range))
end

function Base.checkbounds(seq::BioSequence, locs::AbstractVector{Bool})
    if length(seq) == length(locs)
        return true
    end
    throw(BoundsError(seq, locs))
end

function Base.checkbounds(seq::BioSequence, locs::AbstractVector)
    for i in locs
        checkbounds(seq, i)
    end
    return true
end

function checkdimension(from::Integer, to::Integer)
    if from == to
        return true
    end
    throw(DimensionMismatch(string(
        "attempt to assign ",
        from, " elements to ",
        to,   " elements")))
end

function checkdimension(seq::BioSequence, locs::AbstractVector)
    return checkdimension(length(seq), length(locs))
end

function checkdimension(seq::BioSequence, locs::AbstractVector{Bool})
    return checkdimension(length(seq), sum(locs))
end

# creates a bit mask for given number of bits `n`
mask(n::Integer) = (UInt64(1) << n) - 1
mask(::Type{A}) where {A<:Alphabet} = mask(bitsof(A))

# assumes `i` is positive and `bitsof(A)` is a power of 2
@inline function bitindex(seq::BioSequence{A}, i::Integer) where {A}
    return BitIndex((i + first(seq.part) - 2) << trailing_zeros(bitsof(A)))
end

@inline function inbounds_getindex(seq::BioSequence{A}, i::Integer) where {A}
    j = bitindex(seq, i)
    @inbounds return decode(A, (seq.data[index(j)] >> offset(j)) & mask(A))
end

Base.getindex(seq::BioSequence, part::UnitRange) = BioSequence(seq, part)
Base.view(seq::BioSequence, part::UnitRange) = BioSequence(seq, part)

function Base.setindex!(seq::BioSequence{A},
                        other::BioSequence{A},
                        locs::AbstractVector{<:Integer}) where {A}
    checkbounds(seq, locs)
    checkdimension(other, locs)
    orphan!(seq)
    for (i, x) in zip(locs, other)
        unsafe_setindex!(seq, x, i)
    end
    return seq
end

function Base.setindex!(seq::BioSequence{A},
                        other::BioSequence{A},
                        locs::UnitRange{<:Integer}) where {A}
    checkbounds(seq, locs)
    checkdimension(other, locs)
    return copyto!(seq, locs.start, other, 1)
end

function Base.setindex!(seq::BioSequence{A},
                        other::BioSequence{A},
                        locs::AbstractVector{Bool}) where {A}
    checkbounds(seq, locs)
    checkdimension(other, locs)
    orphan!(seq)
    i = j = 0
    while true
        i = findnext(locs, i + 1)
        if i === nothing
            break
        end
        unsafe_setindex!(seq, other[j+=1], i)
    end
    return seq
end

function Base.setindex!(seq::BioSequence{A}, other::BioSequence{A}, ::Colon) where {A}
    return setindex!(seq, other, 1:lastindex(seq))
end

function Base.setindex!(seq::BioSequence, x, i::Integer)
    checkbounds(seq, i)
    orphan!(seq)
    return unsafe_setindex!(seq, x, i)
end

function Base.setindex!(seq::BioSequence{A}, x, locs::AbstractVector{<:Integer}) where {A}
    checkbounds(seq, locs)
    bin = enc64(seq, x)
    orphan!(seq)
    for i in locs
        encoded_setindex!(seq, bin, i)
    end
    return seq
end

function Base.setindex!(seq::BioSequence{A}, x, locs::AbstractVector{Bool}) where {A}
    checkbounds(seq, locs)
    bin = enc64(seq, x)
    orphan!(seq)
    i = j = 0
    while true
        i = findnext(locs, i + 1)
        if i === nothing
            break
        end
        encoded_setindex!(seq, bin, i)
    end
    return seq
end

Base.setindex!(seq::BioSequence{A}, x, ::Colon) where {A} = setindex!(seq, x, 1:lastindex(seq))

# this is "unsafe" because of no bounds check and no orphan! call
@inline function unsafe_setindex!(seq::BioSequence{A}, x, i::Integer) where {A}
    bin = enc64(seq, x)
    return encoded_setindex!(seq, bin, i)
end

@inline function encoded_setindex!(seq::BioSequence{A}, bin::UInt64, i::Integer) where {A}
    j, r = bitindex(seq, i)
    @inbounds seq.data[j] = (bin << r) | (seq.data[j] & ~(mask(A) << r))
    return seq
end
