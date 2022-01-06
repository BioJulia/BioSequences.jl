###
### LongSequence specific specializations of src/biosequence/indexing.jl
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

# Basic get/setindex methods
@inline function bitindex(x::LongSequence, i::Integer)
    N = BitsPerSymbol(Alphabet(typeof(x)))
    bitindex(N, encoded_data_eltype(typeof(x)), i)
end

firstbitindex(s::SeqOrView) = bitindex(s, firstindex(s))
lastbitindex(s::SeqOrView) = bitindex(s, lastindex(s))

@inline function extract_encoded_element(x::SeqOrView, i::Integer)
    bi = bitindex(x, i % UInt)
    extract_encoded_element(bi, x.data)
end

@inline function encoded_setindex!(s::SeqOrView, v::UInt64, i::BitIndex)
    vi, off = i
    data = s.data
    bits = @inbounds data[vi]
    @inbounds data[vi] = (UInt(v) << off) | (bits & ~(bitmask(i) << off))
    return s
end

# More efficient due to copyto!
function Base.getindex(seq::LongSequence, part::AbstractUnitRange{<:Integer})
	@boundscheck checkbounds(seq, part)
	newseq = typeof(seq)(undef, length(part))
	return copyto!(newseq, 1, seq, first(part), length(part))
end

# More efficient due to copyto!
function Base.setindex!(
    seq::SeqOrView{A},
    other::SeqOrView{A},
    locs::AbstractUnitRange{<:Integer}
) where {A <: Alphabet}
    @boundscheck checkbounds(seq, locs)
    @boundscheck if length(other) != length(locs)
        throw(DimensionMismatch("Attempt to assign $(length(locs)) values to $(length(seq)) destinations"))
    end
    return copyto!(seq, locs.start, other, 1, length(locs))
end

@inline function encoded_setindex!(
    seq::SeqOrView{A},
    bin::Unsigned, 
    i::Integer
) where {A <: Alphabet}
	return encoded_setindex!(seq, bin, bitindex(seq, i))
end

@inline function encoded_setindex!(s::SeqOrView, v::Unsigned, i::BitIndex)
    vi, off = i
    data = s.data
    bits = @inbounds data[vi]
	v_ = v % encoded_data_eltype(typeof(s))
    @inbounds data[vi] = (v_ << off) | (bits & ~(bitmask(i) << off))
    return s
end
