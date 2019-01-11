# Transformations 
# ===============
#
# Methods that manipulate and change a biological sequence.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

"""
    empty!(seq::BioSequence)

Completely empty a biological sequence `seq` of nucleotides.
"""
Base.empty!(seq::BioSequence) = resize!(seq, 0)

"""
    push!(seq::BioSequence, x)

Append a biological symbol `x` to a biological sequence `seq`.
"""
function Base.push!(seq::BioSequence, x)
    resize!(seq, length(seq) + 1)
    unsafe_setindex!(seq, x, lastindex(seq))
    return seq
end

"""
    pop!(seq::BioSequence)

Remove the symbol from the end of a biological sequence `seq` and return it.
Returns a variable of `eltype(seq)`.
"""
function Base.pop!(seq::BioSequence)
    if isempty(seq)
        throw(ArgumentError("sequence must be non-empty"))
    end
    @inbounds x = seq[end]
    deleteat!(seq, lastindex(seq))
    return x
end

"""
    insert!(seq::BioSequence, i, x)

Insert a biological symbol `x` into a biological sequence `seq`, at the given
index `i`.
"""
function Base.insert!(seq::BioSequence, i::Integer, x)
    checkbounds(seq, i)
    resize!(seq, length(seq) + 1)
    copyto!(seq, i + 1, seq, i, lastindex(seq) - i)
    unsafe_setindex!(seq, x, i)
    return seq
end

"""
    deleteat!(seq::BioSequence, range::UnitRange{<:Integer})

Deletes a defined `range` from a biological sequence `seq`.

Modifies the input sequence.
"""
function Base.deleteat!(seq::BioSequence, range::UnitRange{<:Integer})
    checkbounds(seq, range)
    copyto!(seq, range.start, seq, range.stop + 1, length(seq) - range.stop)
    resize!(seq, length(seq) - length(range))
    return seq
end

"""
    deleteat!(seq::BioSequence, i::Integer)

Delete a biological symbol at a single position `i` in a biological sequence
`seq`.

Modifies the input sequence.
"""
function Base.deleteat!(seq::BioSequence, i::Integer)
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
function Base.append!(seq::BioSequence, other::BioSequence)
    resize!(seq, length(seq) + length(other))
    copyto!(seq, lastindex(seq) - length(other) + 1, other, 1)
    return seq
end

"""
    popfirst!(seq)

Remove the symbol from the beginning of a biological sequence `seq` and return
it. Returns a variable of `eltype(seq)`.
"""
function Base.popfirst!(seq::BioSequence)
    if isempty(seq)
        throw(ArgumentError("sequence must be non-empty"))
    end
    @inbounds x = seq[1]
    deleteat!(seq, 1)
    return x
end

"""
    pushfirst!(seq, x)

Insert a biological symbol `x` at the beginning of a biological sequence `seq`.
"""
function Base.pushfirst!(seq::BioSequence, x)
    resize!(seq, length(seq) + 1)
    copyto!(seq, 2, seq, 1, length(seq) - 1)
    unsafe_setindex!(seq, x, 1)
    return seq
end

Base.filter(f::Function, seq::BioSequence) = filter!(f, copy(seq))
Base.map(f::Function, seq::BioSequence) = map!(f, copy(seq))
Base.reverse(seq::BioSequence) = reverse!(copy(seq))

"""
    complement(seq)

Make a complement sequence of `seq`.
"""
function BioSymbols.complement(seq::NucleotideSeq)
    return complement!(copy(seq))
end

"""
    reverse_complement!(seq)

Make a reversed complement sequence of `seq` in place.
"""
function reverse_complement!(seq::NucleotideSeq)
    return complement!(reverse!(seq))
end

"""
    reverse_complement(seq)

Make a reversed complement sequence of `seq`.
"""
function reverse_complement(seq::NucleotideSeq)
    return complement!(reverse(seq))
end

# Shuffle
# -------

function Random.shuffle(seq::BioSequence)
    return shuffle!(copy(seq))
end