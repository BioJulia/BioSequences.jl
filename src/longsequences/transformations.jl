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
    elseif seq.shared
        # May seem wasteful to orphan here, but if you resize! a shared seq,
        # you're probably going to write to it anyway.
        return _orphan!(seq, size)
    else
        if force | (seq_data_len(A, size) > seq_data_len(A, length(seq)))
            resize!(seq.data, seq_data_len(A, size))
        end
        seq.part = 1:size
        return seq
    end
end

function Base.filter!(f::Function, seq::LongSequence{A}) where {A}
    orphan!(seq)

    len = 0
    next = bitindex(seq, 1)
    j = index(next)
    datum::UInt64 = 0
    for i in 1:lastindex(seq)
        x = inbounds_getindex(seq, i)
        if f(x)
            datum |= enc64(seq, x) << offset(next)
            len += 1
            #TODO: Resolve use of bits_per_symbol.
            next += bits_per_symbol(A())
            if index(next) != j
                seq.data[j] = datum
                datum = 0
                j = index(next)
            end
        end
    end
    if offset(next) > 0
        seq.data[j] = datum
    end
    resize!(seq, len)

    return seq
end

function Base.map!(f::Function, seq::LongSequence)
    orphan!(seq)
    for i in 1:lastindex(seq)
        unsafe_setindex!(seq, f(inbounds_getindex(seq, i)), i)
    end
    return seq
end

"""
    reverse!(seq::LongSequence)

Reverse a biological sequence `seq` in place.
"""
Base.reverse!(seq::LongSequence{<:Alphabet}) = _reverse!(orphan!(seq), BitsPerSymbol(seq))

function _reverse!(seq::LongSequence{A}, ::BitsPerSymbol) where {A <: Alphabet}
    i, j = 1, lastindex(seq)
    @inbounds while i < j
        seq[i], seq[j] = seq[j], seq[i]
        i += 1
        j -= 1
    end
    return seq
end

# Inline this large function because it is only called by wrapper functions
@inline function _reverse!(seq::LongSequence{A}, B::BT) where {A <: Alphabet,
    BT <: Union{BitsPerSymbol{2}, BitsPerSymbol{4}, BitsPerSymbol{8}}}

    # Reverse order of chunks and bits in chunks in one pass
    data = seq.data
    len = length(data)
    B = BitsPerSymbol(seq)
    @inbounds for i in 1:len >>> 1
        data[i], data[len-i+1] = reversebits(data[len-i+1], B), reversebits(data[i], B)
    end
    @inbounds if isodd(len)
        data[len >>> 1 + 1] = reversebits(data[len >>> 1 + 1], B)
    end

    # Reversion of chunk bits may have left-shifted data in chunks, so we must
    # shift them back to an offset of zero
    # This is written so it SIMD parallelizes - careful with changes
    lshift = offset(bitindex(seq, last(seq.part)) + bits_per_symbol(A()))
    rshift = 64 - lshift
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

"""
    reverse!(seq::LongSequence)

Reverse a biological sequence.
"""
Base.reverse(seq::LongSequence{A}) where {A<:Alphabet} = _reverse!(copy(seq), BitsPerSymbol(seq))

"""
    complement!(seq)

Make a complement sequence of `seq` in place.
"""
function complement!(seq::LongSequence{A}) where {A<:NucleicAcidAlphabet}
    orphan!(seq)
    next = firstbitindex(seq)
    stop = bitindex(seq, lastindex(seq) + 1)
    seqdata = seq.data
    @inbounds while next < stop
        x = seqdata[index(next)]
        seqdata[index(next)] = complement_bitpar(x, Alphabet(seq))
        next += 64
    end
    return seq
end

###
### Shuffle
###

function Random.shuffle!(seq::LongSequence)
    orphan!(seq) # TODO: Is this call to orphan nessecery, given setindex should call `orphan!` for us?
    # Fisher-Yates shuffle
    for i in 1:lastindex(seq) - 1
        j = rand(i:lastindex(seq))
        seq[i], seq[j] = seq[j], seq[i]
    end
    return seq
end
