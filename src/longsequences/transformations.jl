# Transformations
# ---------------

"""
    push!(seq::LongSequence{A}, x) where {A}

Append a biological symbol `x` to a biological sequence `seq`.
"""
function Base.push!(seq::LongSequence{A}, x) where {A}
    bin = enc64(seq, x)
    resize!(seq, length(seq) + 1)
    encoded_setindex!(seq, bin, lastindex(seq))
    return seq
end

"""
    pop!(seq::LongSequence)

Remove the symbol from the end of a biological sequence `seq` and return it.
Returns a variable of `eltype(seq)`.
"""
function Base.pop!(seq::LongSequence)
    if isempty(seq)
        throw(ArgumentError("sequence must be non-empty"))
    end
    x = seq[end]
    deleteat!(seq, lastindex(seq))
    return x
end

"""
    insert!(seq, i, x)

Insert a biological symbol `x` into a biological sequence `seq`, at the given
index `i`.
"""
function Base.insert!(seq::LongSequence{A}, i::Integer, x) where {A}
    checkbounds(seq, i)
    bin = enc64(seq, x)
    resize!(seq, length(seq) + 1)
    copyto!(seq, i + 1, seq, i, lastindex(seq) - i)
    encoded_setindex!(seq, bin, i)
    return seq
end

"""
    deleteat!(seq::LongSequence, range::UnitRange{<:Integer})

Deletes a defined `range` from a biological sequence `seq`.

Modifies the input sequence.
"""
function Base.deleteat!(seq::LongSequence{A}, range::UnitRange{<:Integer}) where {A}
    checkbounds(seq, range)
    copyto!(seq, range.start, seq, range.stop + 1, length(seq) - range.stop)
    resize!(seq, length(seq) - length(range))
    return seq
end

"""
    deleteat!(seq::LongSequence, i::Integer)

Delete a biological symbol at a single position `i` in a biological sequence
`seq`.

Modifies the input sequence.
"""
function Base.deleteat!(seq::LongSequence, i::Integer)
    checkbounds(seq, i)
    copyto!(seq, i, seq, i + 1, length(seq) - i)
    resize!(seq, length(seq) - 1)
    return seq
end

"""
    append!(seq, other)

Add a biological sequence `other` onto the end of biological sequence `seq`.
Modifies and returns `seq`.
"""
function Base.append!(seq::LongSequence{A},
		      other::LongSequence{A}) where {A}
    resize!(seq, length(seq) + length(other))
    copyto!(seq, lastindex(seq) - length(other) + 1, other, 1)
    return seq
end

"""
    popfirst!(seq)

Remove the symbol from the beginning of a biological sequence `seq` and return
it. Returns a variable of `eltype(seq)`.
"""
function Base.popfirst!(seq::LongSequence)
    if isempty(seq)
        throw(ArgumentError("sequence must be non-empty"))
    end
    x = seq[1]
    deleteat!(seq, 1)
    return x
end

"""
    pushfirst!(seq, x)

Insert a biological symbol `x` at the beginning of a biological sequence `seq`.
"""
function Base.pushfirst!(seq::LongSequence{A}, x) where {A}
    bin = enc64(seq, x)
    resize!(seq, length(seq) + 1)
    copyto!(seq, 2, seq, 1, length(seq) - 1)
    encoded_setindex!(seq, bin, 1)
    return seq
end

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

"""
    empty!(seq)

Completely empty a biological sequence `seq` of nucleotides.
"""
Base.empty!(seq::LongSequence) = resize!(seq, 0)

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

function Base.filter(f::Function, seq::LongSequence)
    return filter!(f, copy(seq))
end

function Base.map!(f::Function, seq::LongSequence)
    orphan!(seq)
    for i in 1:lastindex(seq)
        unsafe_setindex!(seq, f(inbounds_getindex(seq, i)), i)
    end
    return seq
end

function Base.map(f::Function, seq::LongSequence)
    return map!(f, copy(seq))
end

"""
    reverse!(seq)

Reverse a biological sequence `seq` in place.
"""
function Base.reverse!(seq::LongSequence)
    orphan!(seq)
    @inbounds for i in 1:div(lastindex(seq), 2)
	    x = seq[i]
        i′ = lastindex(seq) - i + 1 
	    seq[i] = seq[i′]
	    seq[i′] = x
    end
    return seq
end

"""
    reverse(seq)

Create a sequence which is the reverse of the bioloigcal sequence `seq`.
"""
Base.reverse(seq::LongSequence) = reverse!(copy(seq))

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
function complement!(seq::LongSequence{A}) where {A<:NucleicAcidAlphabet{2}}
    orphan!(seq)
    next = firstbitindex(seq)
    stop = bitindex(seq, lastindex(seq) + 1)
    @inbounds while next < stop
        seq.data[index(next)] = ~seq.data[index(next)]
        next += 64
    end
    return seq
end

"""
    complement!(seq)

Transform `seq` into it's complement.
"""
function complement!(seq::LongSequence{A}) where {A<:NucleicAcidAlphabet{4}}
    orphan!(seq)
    next = firstbitindex(seq)
    stop = bitindex(seq, lastindex(seq) + 1)
    @inbounds while next < stop
        x = seq.data[index(next)]
        seq.data[index(next)] = (
            ((x & 0x1111111111111111) << 3) | ((x & 0x8888888888888888) >> 3) |
            ((x & 0x2222222222222222) << 1) | ((x & 0x4444444444444444) >> 1))
        next += 64
    end
    return seq
end

"""
    complement(seq)

Make a complement sequence of `seq`.
"""
function BioSymbols.complement(seq::LongSequence{A}) where {A<:NucleicAcidAlphabet}
    return complement!(copy(seq))
end

"""
    reverse_complement!(seq)

Make a reversed complement sequence of `seq` in place.
"""
function reverse_complement!(seq::LongSequence{A}) where {A<:NucleicAcidAlphabet}
    return complement!(reverse!(seq))
end

"""
    reverse_complement(seq)

Make a reversed complement sequence of `seq`.
"""
function reverse_complement(seq::LongSequence{A}) where {A<:NucleicAcidAlphabet}
    return complement!(reverse(seq))
end

# Shuffle
# -------

function Random.shuffle(seq::LongSequence)
    return shuffle!(copy(seq))
end

function Random.shuffle!(seq::LongSequence)
    orphan!(seq)
    # Fisher-Yates shuffle
    for i in 1:lastindex(seq)-1
        j = rand(i:lastindex(seq))
        seq[i], seq[j] = seq[j], seq[i]
    end
    return seq
end
