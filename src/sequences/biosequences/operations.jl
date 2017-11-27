# BioSequence operations
# ===================
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

"Remove gap characters from a sequence. Modifies the input sequence."
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


# Predicates
# ----------

"""
    ispalindromic(seq::BioSequence)

Return `true` if `seq` is a palindromic sequence; otherwise return `false`.
"""
function ispalindromic(seq::BioSequence)
    if !(eltype(seq) <: NucleicAcid)
        error("elements must be nucleotide")
    end

    for i in 1:cld(length(seq), 2)
        if seq[i] != complement(seq[end-i+1])
            return false
        end
    end

    return true
end

"""
    hasambiguity(seq::BioSequence)

    Return `true` if `seq` has an ambiguous symbol; otherwise return `false`.
"""
function hasambiguity(seq::BioSequence)
    for x in seq
        if isambiguous(x)
            return true
        end
    end
    return false
end

"""
    isrepetitive(seq::BioSequence, n::Integer=length(seq))

Return `true` if and only if `seq` contains a repetitive subsequence of length `≥ n`.
"""
function isrepetitive(seq::BioSequence, n::Integer=length(seq))
    if n < 0
        error("repetition must be non-negative")
    elseif isempty(seq)
        return n == 0
    end

    rep = 1
    if rep ≥ n
        return true
    end
    last = first(seq)
    for i in 2:endof(seq)
        x = seq[i]
        if x == last
            rep += 1
            if rep ≥ n
                return true
            end
        else
            rep = 1
        end
        last = x
    end

    return false
end
