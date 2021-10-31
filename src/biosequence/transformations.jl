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
    x_ = convert(eltype(seq), x)
    resize!(seq, length(seq) + 1)
    @inbounds seq[end] = x_
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
    @inbounds seq[i] = x
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
    copyto!(seq, lastindex(seq) - length(other) + 1, other, 1, length(other))
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
    @inbounds seq[firstindex(seq)] = x
    return seq
end

Base.filter(f, seq::BioSequence) = filter!(f, copy(seq))

function Base.filter!(f, seq::BioSequence)
    ind = 0
    @inbounds for i in eachindex(seq)
        if f(seq[i])
            ind += 1
        else
            break
        end
    end
    @inbounds for i in ind+1:lastindex(seq)
        v = seq[i]
        if f(v)
            ind += 1
            seq[ind] = v
        end
    end
    return resize!(seq, ind)
end

Base.map(f, seq::BioSequence) = map!(f, copy(seq))

function Base.map!(f, seq::BioSequence)
    @inbounds for i in eachindex(seq)
        seq[i] = f(seq[i])
    end
    seq
end

"""
    reverse(seq::BioSequence)

Create reversed copy of a biological sequence.
"""
Base.reverse(seq::BioSequence) = reverse!(copy(seq))

function Base.reverse!(s::BioSequence)
	i, j = 1, lastindex(s)
	@inbounds while i < j
		s[i], s[j] = s[j], s[i]
		i, j = i + 1, j - 1
	end
	return s
end

"""
    complement(seq)

Make a complement sequence of `seq`.
"""
function BioSymbols.complement(seq::NucleotideSeq)
    return complement!(copy(seq))
end

complement!(seq::NucleotideSeq) = map!(complement, seq)

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

"""
    canonical!(seq::NucleotideSeq)

Transforms the `seq` into its canonical form, if it is not already canonical.
Modifies the input sequence inplace.

For any sequence, there is a reverse complement, which is the same sequence, but
on the complimentary strand of DNA:

```
------->
ATCGATCG
CGATCGAT
<-------
```

!!! note
    Using the [`reverse_complement`](@ref) of a DNA sequence will give give this
    reverse complement.

Of the two sequences, the *canonical* of the two sequences is the lesser of the
two i.e. `canonical_seq < other_seq`.

Using this function on a `seq` will ensure it is the canonical version.
"""
function canonical!(seq::NucleotideSeq)
    if !iscanonical(seq)
        reverse_complement!(seq)
    end
    return seq
end

"""
    canonical(seq::NucleotideSeq)

Create the canonical sequence of `seq`.

"""
canonical(seq::NucleotideSeq) = iscanonical(seq) ? copy(seq) : reverse_complement(seq)

"Create a copy of a sequence with gap characters removed."
ungap(seq::BioSequence)  =  filter(!isgap, seq)

"Remove gap characters from an input sequence."
ungap!(seq::BioSequence) = filter!(!isgap, seq)

###
### Shuffle
###

function Random.shuffle!(seq::BioSequence)
    # Fisher-Yates shuffle
    @inbounds for i in 1:lastindex(seq) - 1
        j = rand(i:lastindex(seq))
        seq[i], seq[j] = seq[j], seq[i]
    end
    return seq
end

function Random.shuffle(seq::BioSequence)
    return shuffle!(copy(seq))
end
