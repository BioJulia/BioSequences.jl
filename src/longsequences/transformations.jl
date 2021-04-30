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
    cp = LongSequence{A}(undef, unsigned(length(seq)))
    reverse_data_copy!(identity, cp.data, seq.data, seq_data_len(seq) % UInt, B)
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
    # We need to account for the fact that the seq may not use all its stored data
    reverse_data!(identity, seq.data, seq_data_len(seq) % UInt, B)
    return zero_offset!(seq)
end


# Reversion of chunk bits may have left-shifted data in chunks, this function right shifts
# all chunks by up to 63 bits.
# This is written so it SIMD parallelizes - careful with changes
@inline function zero_offset!(seq::LongSequence{A}) where A <: Alphabet
    isempty(seq) && return seq
    offs = (64 - offset(bitindex(seq, length(seq)) + bits_per_symbol(A()))) % UInt
    zero_offset!(seq, offs) 
end

@inline function zero_offset!(seq::LongSequence{A}, offs::UInt) where A <: Alphabet
    isempty(seq) && return seq
    rshift = offs
    lshift = 64 - rshift
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
@inline function reverse_data!(pred, data::Vector{UInt64}, len::UInt, B::BT) where {
    BT <: Union{BitsPerSymbol{2}, BitsPerSymbol{4}, BitsPerSymbol{8}}}
    @inbounds @simd ivdep for i in 1:len >>> 1
        data[i], data[len-i+1] = pred(reversebits(data[len-i+1], B)), pred(reversebits(data[i], B))
    end
    @inbounds if isodd(len)
        data[len >>> 1 + 1] = pred(reversebits(data[len >>> 1 + 1], B))
    end
end

@inline function reverse_data_copy!(pred, dst::Vector{UInt64}, src::Vector{UInt64}, len::UInt,
    B::BT) where {BT <: Union{BitsPerSymbol{2}, BitsPerSymbol{4}, BitsPerSymbol{8}}}
    @inbounds @simd for i in eachindex(dst)
        dst[i] = pred(reversebits(src[len - i + 1], B))
    end
end

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

function complement!(s::LongSubSeq{A}) where {A <: NucleicAcidAlphabet}
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
    reverse_data!(pred, seq.data, seq_data_len(seq) % UInt, BitsPerSymbol(seq))
    return zero_offset!(seq)
end

function reverse_complement(seq::LongSequence{<:NucleicAcidAlphabet})
    cp = typeof(seq)(undef, unsigned(length(seq)))
    pred = x -> complement_bitpar(x, Alphabet(seq))
    reverse_data_copy!(pred, cp.data, seq.data, seq_data_len(seq) % UInt, BitsPerSymbol(seq))
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