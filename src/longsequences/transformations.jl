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

"""
    reverse(seq::LongSequence)

Create reversed copy of a biological sequence.
"""
Base.reverse(seq::LongSequence{<:Alphabet}) = _reverse(seq, BitsPerSymbol(seq))

# Fast path for non-inplace reversion
@inline function _reverse(seq::LongSequence{A}, B::BT) where {A <: Alphabet,
    BT <: Union{BitsPerSymbol{2}, BitsPerSymbol{4}, BitsPerSymbol{8}}}
    cp = LongSequence{A}(unsigned(length(seq)))
    reverse_data_copy!(identity, cp.data, seq.data, B)
    return zero_offset!(cp)
end

_reverse(seq::LongSequence{<:Alphabet}, ::BitsPerSymbol) = reverse!(copy(seq))

# Generic fallback
function _reverse!(seq::LongSequence{<:Alphabet}, ::BitsPerSymbol)
    i, j = 1, lastindex(seq)
    @inbounds while i < j
        seq[i], seq[j] = seq[j], seq[i]
        i += 1
        j -= 1
    end
    return seq
end

@inline function _reverse!(seq::LongSequence{<:Alphabet}, B::BT) where {
    BT <: Union{BitsPerSymbol{2}, BitsPerSymbol{4}, BitsPerSymbol{8}}}
    reverse_data!(identity, seq.data, B)
    return zero_offset!(seq)
end

# Reversion of chunk bits may have left-shifted data in chunks, this function right shifts
# all chunks by up to 63 bits.
# This is written so it SIMD parallelizes - careful with changes
@inline function zero_offset!(seq::LongSequence{A}) where A <: Alphabet
    lshift = offset(bitindex(seq, length(seq)) + bits_per_symbol(A()))
    rshift = 64 - lshift
    len = length(seq.data)
    @inbounds if !iszero(lshift)
        this = seq.data[1]
        for i in 1:len-1
            next = seq.data[i+1]
            seq.data[i] = (this >>> (unsigned(rshift) & 63)) | (next << (unsigned(lshift) & 63))
            this = next
        end
        seq.data[len] >>>= (unsigned(rshift) & 63)
    end
    return seq
end

# Reverse chunks in data vector and each symbol within a chunk. Chunks may have nonzero
# offset after use, so use zero_offset!
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

@inline function reverse_data_copy!(pred, dst::Vector{UInt64}, src::Vector{UInt64}, B::BT) where {
    BT <: Union{BitsPerSymbol{2}, BitsPerSymbol{4}, BitsPerSymbol{8}}}
    len = length(dst)
    @inbounds @simd for i in eachindex(dst)
        dst[i] = pred(reversebits(src[len - i + 1], B))
    end
end

"""
    complement!(seq)

Make a complement sequence of `seq` in place.
"""
function complement!(seq::LongSequence{A}) where {A<:NucleicAcidAlphabet}
    next = index(firstbitindex(seq))
    stop = index(lastbitindex(seq))
    seqdata = seq.data
    @inbounds for i in index(firstbitindex(seq)):index(lastbitindex(seq))
        seqdata[i] = complement_bitpar(seqdata[i], Alphabet(seq))
    end
    return seq
end

function reverse_complement!(seq::LongSequence{<:NucleicAcidAlphabet})
    pred = x -> complement_bitpar(x, Alphabet(seq))
    reverse_data!(pred, seq.data, BitsPerSymbol(seq))
    return zero_offset!(seq)
end

function reverse_complement(seq::LongSequence{<:NucleicAcidAlphabet})
    cp = typeof(seq)(unsigned(length(seq)))
    pred = x -> complement_bitpar(x, Alphabet(seq))
    reverse_data_copy!(pred, cp.data, seq.data, BitsPerSymbol(seq))
    return zero_offset!(cp)
end

function Random.shuffle!(seq::LongSequence)
    # Fisher-Yates shuffle
    @inbounds for i in 1:lastindex(seq) - 1
        j = rand(i:lastindex(seq))
        seq[i], seq[j] = seq[j], seq[i]
    end
    return seq
end
