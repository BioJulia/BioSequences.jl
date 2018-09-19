# Indexing
# ========
#
# Indexing methods for mutable biological sequences.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

# assumes `i` is positive and `bitsof(A)` is a power of 2

@inline function bitindex(seq::GeneralSequence, i::Integer)
    return bitindex(BitsPerSymbol(seq), encoded_data_eltype(seq), i + first(seq.part) - 1)
end

@inline function Base.getindex(seq::GeneralSequence, part::UnitRange)
    return GeneralSequence(seq, part)
end
@inline function Base.view(seq::GeneralSequence, part::UnitRange)
    return getindex(seq, part)
end

# Set a single sequence position to a single symbol value.
function Base.setindex!(seq::GeneralSequence, x, i::Integer)
    @boundscheck checkbounds(seq, i)
    orphan!(seq)
    return unsafe_setindex!(seq, x, i)
end

# Set multiple sequence positions to a single symbol value.
function Base.setindex!(seq::GeneralSequence{A}, x,
			locs::AbstractVector{<:Integer}) where {A}

    checkbounds(seq, locs)
    bin = enc64(seq, x)
    orphan!(seq)
    for i in locs
        encoded_setindex!(seq, bin, i)
    end
    return seq
end

function Base.setindex!(seq::GeneralSequence{A}, x,
			locs::AbstractVector{Bool}) where {A}

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

function Base.setindex!(seq::GeneralSequence{A},
                        other::GeneralSequence{A},
                        locs::AbstractVector{<:Integer}) where {A}

    checkbounds(seq, locs)
    checkdimension(other, locs)
    orphan!(seq)
    for (i, x) in zip(locs, other)
        unsafe_setindex!(seq, x, i)
    end
    return seq
end

function Base.setindex!(seq::GeneralSequence{A},
                        other::GeneralSequence{A},
                        locs::UnitRange{<:Integer}) where {A}
    checkbounds(seq, locs)
    checkdimension(other, locs)
    return copyto!(seq, locs.start, other, 1)
end

function Base.setindex!(seq::GeneralSequence{A},
                        other::GeneralSequence{A},
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
        unsafe_setindex!(seq, other[j += 1], i)
    end
    return seq
end

function Base.setindex!(seq::GeneralSequence{A},
			other::GeneralSequence{A}, ::Colon) where {A}
    return setindex!(seq, other, 1:lastindex(seq))
end

function Base.setindex!(seq::GeneralSequence{A}, x, ::Colon) where {A}
    return setindex!(seq, x, 1:lastindex(seq))
end

# this is "unsafe" because of no bounds check and no orphan! call
@inline function unsafe_setindex!(seq::GeneralSequence{A}, x, i::Integer) where {A}
    bin = enc64(seq, x)
    return encoded_setindex!(seq, bin, i)
end

function enc64(::GeneralSequence{A}, x) where {A}
    #TODO: Resolve these two use cases of A().
    return UInt64(encode(A(), convert(eltype(A()), x)))
end

@inline function encoded_setindex!(seq::GeneralSequence{A},
				   bin::UInt64, i::Integer) where {A}
    j, r = bitindex(seq, i)
    data = encoded_data(seq)
    @inbounds data[j] = (bin << r) | (data[j] & ~(bindata_mask(seq) << r))
    return seq
end
