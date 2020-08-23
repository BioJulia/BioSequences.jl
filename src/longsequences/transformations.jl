###
### LongSequence specific specializations of src/biosequence/transformations.jl
###

"""
    resize!(seq, size, [force::Bool])

Resize a biological sequence `seq`, to a given `size`. Does not resize the underlying data
array unless the new size does not fit. If `force`, always resize underlying data array.
"""
function Base.resize!(seq::LongSequence{A}, size::Integer, force::Bool=false) where {A}
    if size < 0
        throw(ArgumentError("size must be non-negative"))
    else
        if force | (seq_data_len(A, size) > seq_data_len(A, length(seq)))
            resize!(seq.data, seq_data_len(A, size))
        end
        seq.len = size
        return seq
    end
end

function Base.filter!(f, seq::LongSequence)
    writeindex = bitindex(seq, 0)
    chunk = zero(UInt64)
    bps = bits_per_symbol(seq)
    for i in eachindex(seq)
        encoded_symbol = extract_encoded_element(bitindex(seq, i), encoded_data(seq))
        symbol = decode(Alphabet(seq), encoded_symbol)
        if f(symbol)
            writeindex += bps
            chunk |= encoded_symbol << offset(writeindex)
            # If the offset is now zero, the chunk is full
            if iszero(offset(writeindex + bps))
                seq.data[index(writeindex)] = chunk
                chunk = zero(UInt64)
            end
        end
    end
    # If it's not zero, we need to write the first chunk
    if !iszero(offset(writeindex + bps))
        seq.data[index(writeindex)] = chunk
    end
    resize!(seq, div((writeindex + bps).val, bps))
end

function Base.map!(f, seq::LongSequence)
    for i in 1:lastindex(seq)
        unsafe_setindex!(seq, f(inbounds_getindex(seq, i)), i)
    end
    return seq
end

"""
    reverse!(seq::LongSequence)

Reverse a biological sequence `seq` in place.
"""
Base.reverse!(seq::LongSequence{<:Alphabet}) = _reverse!(seq, BitsPerSymbol(seq))

# Specialized methods for 2, 4, 8 bits per symbol
@inline function _reverse!(seq::LongSequence{<:Alphabet}, B::BT) where {
    BT <: Union{BitsPerSymbol{2}, BitsPerSymbol{4}, BitsPerSymbol{8}}}
    reverse_data!(identity, seq.data, B)
	rightshift!(seq.data, reverse_offset(seq))
    return seq
end

_reverse!(seq::LongSequence, ::BitsPerSymbol) = invoke(reverse!, Tuple{BioSequence}, seq)

"""
    complement!(seq)

Make a complement sequence of `seq` in place.
"""
function complement!(seq::LongSequence{A}) where {A<:NucleicAcidAlphabet}
    seqdata = seq.data
    @inbounds for i in eachindex(seqdata)
        seqdata[i] = complement_bitpar(seqdata[i], Alphabet(seq))
    end
    return seq
end

function complement!(s::SeqView{A}) where {A <: NucleicAcidAlphabet}
	bps = bits_per_symbol(A())
	bi = firstbitindex(s)
	i = 1
	stop = lastbitindex(s) + bps
	@inbounds while (!iszero(offset(bi)) & (bi < stop))
		s[i] = complement(s[i])
		bi += bps
		i += 1
	end
	@inbounds for j in index(bi):index(stop)-1
		s.data[j] = complement_bitpar(s.data[j], Alphabet(s))
		bi += 64
		i += symbols_per_data_element(s)
	end
	@inbounds while bi < stop
		s[i] = complement(s[i])
		bi += bps
		i += 1
	end
	return s
end

function reverse_complement!(seq::LongSequence{<:NucleicAcidAlphabet})
    pred = x -> complement_bitpar(x, Alphabet(seq))
    reverse_data!(pred, seq.data, BitsPerSymbol(seq))
	rightshift!(seq.data, reverse_offset(seq))
    return seq
end

function reverse_complement(seq::LongSequence{<:NucleicAcidAlphabet})
	return reverse_complement!(copy(seq))
end

#-- utility functions

function reverse_offset(seq::LongSequence)
	return (64 - offset(bitindex(seq, lastindex(seq) + 1))) & 63
end

# This is written so it SIMD parallelizes - careful with changes
@inline function rightshift!(v::Vector{UInt64}, offset::Integer)
	len = length(v)
    @inbounds if !iszero(offset)
        this = v[1]
        for i in 1:len-1
            next = v[i+1]
            v[i] = (this >>> (offset & 63)) | (next << ((64-offset) & 63))
            this = next
        end
        v[len] >>>= (offset & 63)
    end
    return v
end

# Reverse chunks in data vector and each symbol within a chunk. Chunks may have nonzero
# offset after use, so remember to use rightshift!
@inline function reverse_data!(pred, data::Vector{UInt64}, B::BT) where {
    BT <: Union{BitsPerSymbol{2}, BitsPerSymbol{4}, BitsPerSymbol{8}}}
    len = length(data)
    @inbounds @simd ivdep for i in 1:len >>> 1
        data[i], data[len-i+1] = pred(reversebits(data[len-i+1], B)), pred(reversebits(data[i], B))
    end
    @inbounds if isodd(len)
        data[len >>> 1 + 1] = pred(reversebits(data[len >>> 1 + 1], B))
    end
end
