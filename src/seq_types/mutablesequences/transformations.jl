# Transformations
# ---------------

"""
    push!(seq::MutableBioSequence{A}, x) where {A}

Append a biological symbol `x` to a biological sequence `seq`.
"""
function Base.push!(seq::MutableBioSequence{A}, x) where {A}
    bin = enc64(seq, x)
    resize!(seq, length(seq) + 1)
    encoded_setindex!(seq, bin, lastindex(seq))
    return seq
end

"""
    pop!(seq::MutableBioSequence)

Remove the symbol from the end of a biological sequence `seq` and return it.
Returns a variable of `eltype(seq)`.
"""
function Base.pop!(seq::MutableBioSequence)
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
function Base.insert!(seq::MutableBioSequence{A}, i::Integer, x) where {A}
    checkbounds(seq, i)
    bin = enc64(seq, x)
    resize!(seq, length(seq) + 1)
    copyto!(seq, i + 1, seq, i, lastindex(seq) - i)
    encoded_setindex!(seq, bin, i)
    return seq
end

"""
    deleteat!(seq::MutableBioSequence, range::UnitRange{<:Integer})

Deletes a defined `range` from a biological sequence `seq`.

Modifies the input sequence.
"""
function Base.deleteat!(seq::MutableBioSequence{A}, range::UnitRange{<:Integer}) where {A}
    checkbounds(seq, range)
    copyto!(seq, range.start, seq, range.stop + 1, length(seq) - range.stop)
    resize!(seq, length(seq) - length(range))
    return seq
end

"""
    deleteat!(seq::MutableBioSequence, i::Integer)

Delete a biological symbol at a single position `i` in a biological sequence
`seq`.

Modifies the input sequence.
"""
function Base.deleteat!(seq::MutableBioSequence, i::Integer)
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
function Base.append!(seq::MutableBioSequence{A},
		      other::MutableBioSequence{A}) where {A}
    resize!(seq, length(seq) + length(other))
    copyto!(seq, lastindex(seq) - length(other) + 1, other, 1)
    return seq
end

"""
    popfirst!(seq)

Remove the symbol from the beginning of a biological sequence `seq` and return
it. Returns a variable of `eltype(seq)`.
"""
function Base.popfirst!(seq::MutableBioSequence)
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
function Base.pushfirst!(seq::MutableBioSequence{A}, x) where {A}
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
function Base.resize!(seq::MutableBioSequence{A}, size::Integer) where {A}
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
Base.empty!(seq::MutableBioSequence) = resize!(seq, 0)

function Base.filter!(f::Function, seq::MutableBioSequence{A}) where {A}
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
            next += bits_per_symbol(A)
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

function Base.filter(f::Function, seq::MutableBioSequence)
    return filter!(f, copy(seq))
end

function Base.map!(f::Function, seq::MutableBioSequence)
    orphan!(seq)
    for i in 1:lastindex(seq)
        unsafe_setindex!(seq, f(inbounds_getindex(seq, i)), i)
    end
    return seq
end

function Base.map(f::Function, seq::MutableBioSequence)
    return map!(f, copy(seq))
end

"""
    reverse!(seq)

Reverse a biological sequence `seq` in place.
"""
function Base.reverse!(seq::MutableBioSequence)
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
Base.reverse(seq::MutableBioSequence) = reverse!(copy(seq))

@generated function Base.reverse(seq::MutableBioSequence{A}) where {A<:NucleicAcidAlphabet}
    n = bits_per_symbol(A)
    if n == 2
        nucrev = :nucrev2
    elseif n == 4
        nucrev = :nucrev4
else
        error("n (= $n) ∉ (2, 4)")
    end

    quote
        data = Vector{UInt64}(undef, seq_data_len(A, length(seq)))
        i = 1
        next = BitIndex(seq, lastindex(seq))
        stop = BitIndex(seq, 0)
        r = rem(offset(next) + $n, 64)
        if r == 0
            @inbounds while next - stop > 0
                x = seq.data[index(next)]
                data[i] = $nucrev(x)
                i += 1
                next -= 64
            end
        else
            @inbounds while next - stop > 64
                j = index(next)
                x = (seq.data[j] << (64 - r)) | (seq.data[j-1] >> r)
                data[i] = $nucrev(x)
                i += 1
                next -= 64
            end
            if next - stop > 0
                j = index(next)
                x = seq.data[j] << (64 - r)
                if r < next - stop
                    x |= seq.data[j-1] >> r
                end
                data[i] = $nucrev(x)
            end
        end
        return MutableBioSequence{A}(data, 1:length(seq), false)
    end
end

@inline function nucrev2(x::UInt64)
     x = (x & 0x3333333333333333) <<  2 | (x & 0xCCCCCCCCCCCCCCCC) >>  2
     x = (x & 0x0F0F0F0F0F0F0F0F) <<  4 | (x & 0xF0F0F0F0F0F0F0F0) >>  4
     x = (x & 0x00FF00FF00FF00FF) <<  8 | (x & 0xFF00FF00FF00FF00) >>  8
     x = (x & 0x0000FFFF0000FFFF) << 16 | (x & 0xFFFF0000FFFF0000) >> 16
     x = (x & 0x00000000FFFFFFFF) << 32 | (x & 0xFFFFFFFF00000000) >> 32
     return x
end

@inline function nucrev4(x::UInt64)
     x = (x & 0x0F0F0F0F0F0F0F0F) <<  4 | (x & 0xF0F0F0F0F0F0F0F0) >>  4
     x = (x & 0x00FF00FF00FF00FF) <<  8 | (x & 0xFF00FF00FF00FF00) >>  8
     x = (x & 0x0000FFFF0000FFFF) << 16 | (x & 0xFFFF0000FFFF0000) >> 16
     x = (x & 0x00000000FFFFFFFF) << 32 | (x & 0xFFFFFFFF00000000) >> 32
     return x
end

"""
    complement!(seq)

Make a complement sequence of `seq` in place.
"""
function complement!(seq::MutableBioSequence{A}) where {A<:NucleicAcidAlphabet{2}}
    orphan!(seq)
<<<<<<< HEAD
    next = BitIndex(seq, 1)
    stop = BitIndex(seq, lastindex(seq) + 1)
=======
    next = bitindex(seq, 1)
    stop = bitindex(seq, endof(seq) + 1)
>>>>>>> Parametric BitIndex (#4)
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
function complement!(seq::MutableBioSequence{A}) where {A<:NucleicAcidAlphabet{4}}
    orphan!(seq)
<<<<<<< HEAD
    next = BitIndex(seq, 1)
    stop = BitIndex(seq, lastindex(seq) + 1)
=======
    next = bitindex(seq, 1)
    stop = bitindex(seq, endof(seq) + 1)
>>>>>>> Parametric BitIndex (#4)
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
function BioSymbols.complement(seq::MutableBioSequence{A}) where {A<:NucleicAcidAlphabet}
    return complement!(copy(seq))
end

"""
    reverse_complement!(seq)

Make a reversed complement sequence of `seq` in place.
"""
function reverse_complement!(seq::MutableBioSequence{A}) where {A<:NucleicAcidAlphabet}
    return complement!(reverse!(seq))
end

"""
    reverse_complement(seq)

Make a reversed complement sequence of `seq`.
"""
function reverse_complement(seq::MutableBioSequence{A}) where {A<:NucleicAcidAlphabet}
    return complement!(reverse(seq))
end

# Shuffle
# -------

function Random.shuffle(seq::MutableBioSequence)
    return shuffle!(copy(seq))
end

function Random.shuffle!(seq::MutableBioSequence)
    orphan!(seq)
    # Fisher-Yates shuffle
    for i in 1:lastindex(seq)-1
        j = rand(i:lastindex(seq))
        seq[i], seq[j] = seq[j], seq[i]
    end
    return seq
end
