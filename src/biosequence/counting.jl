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

function count_naive(f::Function, seq::BioSequence, args...)
    n = 0
    @inbounds for x in seq
        if f(x, args...)
            n += 1
        end
    end
    return n
end

"""
Count how many positions in a sequence satisfy a condition (i.e. f(seq[i]) -> true).

The first argument should be a function which accepts an element of the sequence
as its first parameter, additional arguments may be passed with `args...`.
"""
Base.count(f::Function, seq::BioSequence, args...) = count_naive(f, seq, args...)

###
### GC content
###

"""
    gc_content(seq::BioSequence)

Calculate GC content of `seq`.
"""
gc_content(seq::NucleotideSeq) = isempty(seq) ? 0.0 : count_gc(seq) / length(seq)

function count_gc(seq::BioSequence)
    return count(isGC, seq)
end