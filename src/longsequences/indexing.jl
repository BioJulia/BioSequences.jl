###
### LongSequence specific specializations of src/biosequence/indexing.jl
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

# assumes `i` is positive and `bitsof(A)` is a power of 2

@inline function bitindex(seq::LongSequence, i::Integer)
    return bitindex(BitsPerSymbol(seq), encoded_data_eltype(seq), i + first(seq.part) - 1)
end

@inline function Base.getindex(seq::LongSequence, part::UnitRange)
    return LongSequence(seq, part)
end
@inline function Base.view(seq::LongSequence, part::UnitRange)
    return getindex(seq, part)
end

# Set a single sequence position to a single symbol value.
function Base.setindex!(seq::LongSequence, x, i::Integer)
    @boundscheck checkbounds(seq, i)
    orphan!(seq)
    return unsafe_setindex!(seq, x, i)
end

# Set multiple sequence positions to a single symbol value.
function Base.setindex!(seq::LongSequence{A}, x, locs::AbstractVector{<:Integer}) where {A}
    @boundscheck checkbounds(seq, locs)
    orphan!(seq)
    return unsafe_setindex!(seq, x, locs)
end

function Base.setindex!(seq::LongSequence{A}, x, locs::AbstractVector{Bool}) where {A}
    @boundscheck checkbounds(seq, locs)
    orphan!(seq)
    return unsafe_setindex!(seq, x, locs)
end

function Base.setindex!(seq::LongSequence{A},
                        other::LongSequence{A},
                        locs::AbstractVector{<:Integer}) where {A}
    @boundscheck checkbounds(seq, locs)
    checkdimension(other, locs)
    orphan!(seq)
    return unsafe_setindex!(seq, other, locs)
end

function Base.setindex!(seq::LongSequence{A},
                        other::LongSequence{A},
                        locs::AbstractVector{Bool}) where {A}
    @boundscheck checkbounds(seq, locs)
    checkdimension(other, locs)
    orphan!(seq)
    return unsafe_setindex!(seq, other, locs)
end

function Base.setindex!(seq::LongSequence{A},
                        other::LongSequence{A},
                        locs::UnitRange{<:Integer}) where {A}
    @boundscheck checkbounds(seq, locs)
    checkdimension(other, locs)
    return copyto!(seq, locs.start, other, 1, length(locs))
end

function Base.setindex!(seq::LongSequence{A},
			            other::LongSequence{A}, ::Colon) where {A}
    return setindex!(seq, other, 1:lastindex(seq))
end

function Base.setindex!(seq::LongSequence{A}, x, ::Colon) where {A}
    return setindex!(seq, x, 1:lastindex(seq))
end

# These are "unsafe" because of no bounds check and no orphan! call
@inline function unsafe_setindex!(seq::LongSequence{A}, x, i::Integer) where {A}
    bin = enc64(seq, x)
    return encoded_setindex!(seq, bin, i)
end

@inline function unsafe_setindex!(seq::LongSequence{A}, x, locs::AbstractVector{<:Integer}) where {A}
    bin = enc64(seq, x)
    for i in locs
        encoded_setindex!(seq, bin, i)
    end
    return seq
end

function unsafe_setindex!(seq::LongSequence{A}, x, locs::AbstractVector{Bool}) where {A}
    bin = enc64(seq, x)
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

function unsafe_setindex!(seq::LongSequence{A}, other::LongSequence{A}, locs::AbstractVector{<:Integer}) where{A}
    for (i, x) in zip(locs, other)
        unsafe_setindex!(seq, x, i)
    end
    return seq
end

function unsafe_setindex!(seq::LongSequence{A},
                          other::LongSequence{A},
                          locs::AbstractVector{Bool}) where {A}
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

function enc64(::LongSequence{A}, x) where {A}
    #TODO: Resolve these two use cases of A().
    return UInt64(encode(A(), convert(eltype(A()), x)))
end

@inline function encoded_setindex!(seq::LongSequence{A},
				   bin::UInt64, i::Integer) where {A}
    j, r = bitindex(seq, i)
    data = encoded_data(seq)
    @inbounds data[j] = (bin << r) | (data[j] & ~(bindata_mask(seq) << r))
    return seq
end

#=
function Base.iterate(seq::LongSequence{A}, off::Int=first(seq.part)-1) where {A <: Alphabet}
    off == last(seq.part) && return nothing
    bps = bits_per_symbol(A())
    @inbounds chunk = seq.data[(off >>> index_shift(BitsPerSymbol(A()))) + 1]
    shift = (unsigned(off) * bps) & 63
    encoding = (chunk >>> shift) & bitmask(A())
    return decode(A(), encoding), off + 1
end

=#
