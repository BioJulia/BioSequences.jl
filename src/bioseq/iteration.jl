# Iteration
# =========
#
# Types and methods for iterating over biological sequences.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

# Ambiguous nucleotides iterator
# ------------------------------

struct AmbiguousNucleicAcidIterator{A<:NucAlphs}
    seq::BioSequence{A}
end

ambiguous_positions(seq::BioSequence) = AmbiguousNucleicAcidIterator(seq)

Base.iterate(it::AmbiguousNucleicAcidIterator) = iterate(it, find_next_ambiguous(it.seq, 1))

function Base.iterate(it::AmbiguousNucleicAcidIterator, nextpos::Int)
    if nextpos == 0
        return nothing
    else
        nextpos, find_next_ambiguous(it.seq, nextpos + 1)
    end
end

Base.IteratorSize(::AmbiguousNucleicAcidIterator) = Base.SizeUnknown()

function find_next_ambiguous(seq::BioSequence{A}, i::Integer) where {A<:TwoBitNucs}
    # no ambiguity
    return 0
end

function find_next_ambiguous(seq::BioSequence{A}, from::Integer) where {A<:FourBitNucs}
    for i in max(from, 1):lastindex(seq)
        nt = inbounds_getindex(seq, i)
        if isambiguous(nt)
            return i
        end
    end
    # no ambiguity
    return 0
end
