###
### LongSequence specific specializations of src/biosequence/indexing.jl
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

# Internally, LongSequence and LongSubSeq are indexed by a bitindex
function bitindex(x::LongSequence, i::Integer)
    N = BitsPerSymbol(Alphabet(typeof(x)))
    bitindex(N, encoded_data_eltype(typeof(x)), i)
end

function bitindex(x::LongSubSeq, i::Integer)
    N = BitsPerSymbol(Alphabet(typeof(x)))
    bitindex(N, encoded_data_eltype(typeof(x)), i - first(x.part) + 1)
end

extract_encoded_element(x::SeqOrView, i::Integer) = extract_encoded_element(x, bitindex(x, i))

function extract_encoded_element(x::SeqOrView, i::BitIndex{N, UInt64}) where N
    extract_encoded_element(i, x.data)
end

####################

# More efficient due to copyto!
function Base.getindex(seq::LongSequence, part::UnitRange{<:Integer})
	@boundscheck checkbounds(seq, part)
	newseq = typeof(seq)(undef, length(part))
	return copyto!(newseq, 1, seq, first(part), length(part))
end

# More efficient due to copyto!
function Base.setindex!(seq::SeqOrView{A},
                        other::SeqOrView{A},
                        locs::UnitRange{<:Integer}) where {A <: Alphabet}
    @boundscheck checkbounds(seq, locs)
    checkdimension(other, locs)
    return copyto!(seq, locs.start, other, 1, length(locs))
end

@inline function encoded_setindex!(seq::SeqOrView{A},
				   bin::Unsigned, i::Integer) where {A <: Alphabet}
	return encoded_setindex!(seq, bin, bitindex(seq, i))
end

@inline function encoded_setindex!(s::SeqOrView, v::Unsigned, i::BitIndex)
    vi, off = i
    data = encoded_data(s)
    bits = @inbounds data[vi]
	v_ = v % encoded_data_eltype(s)
    @inbounds data[vi] = (v_ << off) | (bits & ~(bindata_mask(s) << off))
    return s
end

function Base.setindex!(seq::SeqOrView, x, locs::AbstractVector{<:Integer})
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
                        other::SeqOrView{A},
                        locs::AbstractVector{<:Integer}) where {A <: Alphabet}
    @boundscheck checkbounds(seq, locs)
    checkdimension(other, locs)
	@inbounds for (i, n) in zip(locs, other)
		seq[i] = n
	end
	return seq
end
