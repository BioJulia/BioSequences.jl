# Indexing
# ========
#
# Indexing methods for mutable biological sequences.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

# assumes `i` is positive and `bitsof(A)` is a power of 2
@inline function BitIndex(seq::MutableBioSequence, i::Integer)
    nbits = bitsof(alphabet_t(seq))
    return BitIndex((i + first(seq.part) - 2) << trailing_zeros(nbits))
end

@inline function Base.getindex(seq::MutableBioSequence, part::UnitRange)
    return MutableBioSequence(seq, part)
end
@inline function Base.view(seq::MutableBioSequence, part::UnitRange)
    return getindex(seq, part)
end

# Set a single sequence position to a single symbol value.
function Base.setindex!(seq::MutableBioSequence, x, i::Integer)
    @boundscheck checkbounds(seq, i)
    orphan!(seq)
    return unsafe_setindex!(seq, x, i)
end

# Set multiple sequence positions to a single symbol value.
function Base.setindex!(seq::MutableBioSequence{A}, x, 
			locs::AbstractVector{<:Integer}) where {A}

    checkbounds(seq, locs)
    bin = enc64(seq, x)
    orphan!(seq)
    for i in locs
        encoded_setindex!(seq, bin, i)
    end
    return seq
end

function Base.setindex!(seq::MutableBioSequence{A}, x,
			locs::AbstractVector{Bool}) where {A}

    checkbounds(seq, locs)
    bin = enc64(seq, x)
    orphan!(seq)
    i = j = 0
    while (i = findnext(locs, i + 1)) > 0
        encoded_setindex!(seq, bin, i)
    end
    return seq
end

function Base.setindex!(seq::MutableBioSequence{A},
                        other::MutableBioSequence{A},
                        locs::AbstractVector{<:Integer}) where {A}
    
    checkbounds(seq, locs)
    checkdimension(other, locs)
    orphan!(seq)
    for (i, x) in zip(locs, other)
        unsafe_setindex!(seq, x, i)
    end
    return seq
end

function Base.setindex!(seq::MutableBioSequence{A},
                        other::MutableBioSequence{A},
                        locs::UnitRange{<:Integer}) where {A}
    checkbounds(seq, locs)
    checkdimension(other, locs)
    return copy!(seq, locs.start, other, 1)
end

function Base.setindex!(seq::MutableBioSequence{A},
                        other::MutableBioSequence{A},
                        locs::AbstractVector{Bool}) where {A}
    checkbounds(seq, locs)
    checkdimension(other, locs)
    orphan!(seq)
    i = j = 0
    while (i = findnext(locs, i + 1)) > 0
        unsafe_setindex!(seq, other[j += 1], i)
    end
    return seq
end

function Base.setindex!(seq::MutableBioSequence{A},
			other::MutableBioSequence{A}, ::Colon) where {A}
    return setindex!(seq, other, 1:endof(seq))
end

function Base.setindex!(seq::MutableBioSequence{A}, x, ::Colon) where {A}
    return setindex!(seq, x, 1:endof(seq))
end

# this is "unsafe" because of no bounds check and no orphan! call
@inline function unsafe_setindex!(seq::MutableBioSequence{A}, x, i::Integer) where {A}
    bin = enc64(seq, x)
    return encoded_setindex!(seq, bin, i)
end

function enc64(::MutableBioSequence{A}, x) where {A}
    return UInt64(encode(A, convert(eltype(A), x)))
end

@inline function encoded_setindex!(seq::MutableBioSequence{A},
				   bin::UInt64, i::Integer) where {A}
    j, r = BitIndex(seq, i)
    data = bindata(seq)
    @inbounds data[j] = (bin << r) | (data[j] & ~(bindata_mask(seq) << r))
    return seq
end
