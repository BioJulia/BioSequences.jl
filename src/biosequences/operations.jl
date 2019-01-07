# BioSequence operations
# ======================
#
# Generalised operations on biological sequence types.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

"""
Count how many nucleotides satisfy a condition (i.e. f(seq[i]) -> true).

The first argument should be a function which accepts a nucleotide as its parameter.
"""
function Base.count(f::Function, seq::BioSequence)
    n = 0
    @inbounds for x in seq
        if f(x)
            n += 1
        end
    end
    return n
end

# Transformations
# ---------------

"Create a copy of a sequence with gap characters removed."
ungap(seq::BioSequence)  =  filter(x -> x != gap(eltype(seq)), seq)

"Remove gap characters from an input sequence."
ungap!(seq::BioSequence) = filter!(x -> x != gap(eltype(seq)), seq)


# GC content
# ----------

"""
    gc_content(seq::BioSequence)

Calculate GC content of `seq`.
"""
function gc_content(seq::BioSequence)
    if !(eltype(seq) <: NucleicAcid)
        throw(ArgumentError("Not a nucleic acid sequence"))
    end
    if isempty(seq)
        return 0.0
    else
        return count_gc(seq) / length(seq)
    end
end

function count_gc(seq::BioSequence)
    return count(isGC, seq)
end
