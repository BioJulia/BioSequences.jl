# Transformations
# ---------------

"""
    resize!(seq, size)

Resize a biological sequence `seq`, to a given `size`.
"""
function Base.resize!(seq::LongSequence{A}, size::Integer) where {A}
    if size < 0
        throw(ArgumentError("size must be non-negative"))
    end
    orphan!(seq, size)
    resize!(seq.data, seq_data_len(A, size + seq.part.start - 1))
    seq.part = seq.part.start:seq.part.start+size-1
    return seq
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
function Base.reverse!(seq::LongSequence)
    orphan!(seq) # TODO: Is the orphan call really nessecery given the indexing calls will also call orphan?
    @inbounds for i in 1:div(lastindex(seq), 2)
	    x = seq[i]
        i′ = lastindex(seq) - i + 1 
	    seq[i] = seq[i′]
	    seq[i′] = x
    end
    return seq
end

function Base.reverse(seq::LongSequence{A}) where {A<:NucleicAcidAlphabet}
    data = Vector{UInt64}(undef, seq_data_len(A, length(seq)))
    i = 1
    next = lastbitindex(seq)
    stop = bitindex(seq, 0)
    r = rem(offset(next) + bits_per_symbol(seq), 64)
    if r == 0
        @inbounds while next - stop > 0
            x = seq.data[index(next)]
            data[i] = reversebits(x, BitsPerSymbol(seq))
            i += 1
            next -= 64
        end
    else
        @inbounds while next - stop > 64
            j = index(next)
            x = (seq.data[j] << (64 - r)) | (seq.data[j - 1] >> r)
            data[i] = reversebits(x, BitsPerSymbol(seq))
            i += 1
            next -= 64
        end
        if next - stop > 0
            j = index(next)
            x = seq.data[j] << (64 - r)
            if r < next - stop
                x |= seq.data[j - 1] >> r
            end
            data[i] = reversebits(x, BitsPerSymbol(seq))
        end
    end
    return LongSequence{A}(data, 1:length(seq), false)
end

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
