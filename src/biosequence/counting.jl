###
### Counting
###
### Counting operations on biological sequence types.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

###
### Naive counting
###

function count_naive(pred, seq::BioSequence)
    n = 0
    @inbounds for i in eachindex(seq)
        n += pred(seq[i])::Bool
    end
    return n
end

function count_naive(pred, seqa::BioSequence, seqb::BioSequence)
    n = 0
    @inbounds for i in 1:min(length(seqa), length(seqb))
        n += pred(seqa[i], seqb[i])::Bool
    end
    return n
end

"""
Count how many positions in a sequence satisfy a condition (i.e. f(seq[i]) -> true).

The first argument should be a function which accepts an element of the sequence
as its first parameter, additional arguments may be passed with `args...`.
"""
Base.count(pred, seq::BioSequence) = count_naive(pred, seq)

Base.count(pred, seqa::BioSequence, seqb::BioSequence) = count_naive(pred, seqa, seqb)

###
### Aliases for various uses of `count`.
###

"""
    gc_content(seq::BioSequence)

Calculate GC content of `seq`.
"""
gc_content(seq::NucleotideSeq) = isempty(seq) ? 0.0 : count(isGC, seq) / length(seq)

mismatches(seqa::BioSequence, seqb::BioSequence) = count(!=, seqa, seqb)
matches(seqa::BioSequence, seqb::BioSequence) = count(==, seqa, seqb)