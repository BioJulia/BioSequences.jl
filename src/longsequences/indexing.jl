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

function Base.setindex!(seq::SeqOrView, x, locs::AbstractVector{<:Integer})
	@boundscheck checkbounds(seq, locs)
	unsafe_setindex!(seq, x, locs)
end

function Base.setindex!(seq::SeqOrView, x::BioSequence, locs::AbstractVector{<:Integer})
	@boundscheck checkbounds(seq, locs)
	unsafe_setindex!(seq, x, locs)
end

@inline function unsafe_setindex!(seq::SeqOrView, x, locs::AbstractVector{<:Integer})
	bin = encode(Alphabet(seq), convert(eltype(seq), x))
    for i in locs
        encoded_setindex!(seq, bin, bitindex(seq, i))
    end
    return seq
end

# Only to avoid ambiguity
function Base.setindex!(seq::SeqOrView, x::SeqOrView, locs::AbstractVector{Bool})
	@boundscheck checkbounds(seq, locs)
	checkdimension(x, locs)
	j = 0
	@inbounds for i in eachindex(locs)
		if locs[i]
			j += 1
			seq[i] = x[j]
		end
	end
	return seq
end

# To save encoding.
function unsafe_setindex!(seq::SeqOrView, x, locs::AbstractVector{Bool})
    bin = encode(Alphabet(seq), convert(eltype(seq), x))
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

# To avoid ambiguity errors
function Base.setindex!(seq::SeqOrView{A},
                        other::BioSequence{A},
                        locs::AbstractVector{<:Integer}) where {A <: Alphabet}
    @boundscheck checkbounds(seq, locs)
    checkdimension(other, locs)
    unsafe_setindex!(seq, other, locs)
end

function unsafe_setindex!(seq::SeqOrView{A},
                        other::BioSequence,
                        locs::AbstractVector{<:Integer}) where {A <: Alphabet}
	@inbounds for (i, n) in zip(locs, other)
		seq[i] = n
	end
	return seq
end
