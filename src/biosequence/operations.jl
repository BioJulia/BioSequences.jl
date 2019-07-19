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

###
### Finders
###

function Base.findnext(val, seq::BioSequence, start::Integer)
    checkbounds(seq, start)
    v = convert(eltype(seq), val)
    for i in start:lastindex(seq)
        if inbounds_getindex(seq, i) == v
            return i
        end
    end
    return nothing
end

function Base.findnext(f::Function, seq::BioSequence, start::Integer)
    checkbounds(seq, start)
    for i in start:lastindex(seq)
        if f(inbounds_getindex(seq, i))
            return i
        end
    end
    return nothing
end

# No ambiguous sites can exist in a nucleic acid sequence using the two-bit alphabet.
Base.findnext(f::typeof(isambiguous), seq::BioSequence{<:NucleicAcidAlphabet{2}}, from::Integer) = nothing

function Base.findprev(val, seq::BioSequence, start::Integer)
    checkbounds(seq, start)
    v = convert(eltype(seq), val)
    for i in start:-1:1
        if inbounds_getindex(seq, i) == v
            return i
        end
    end
    return nothing
end

function Base.findprev(f::Function, seq::BioSequence, start::Integer)
    checkbounds(seq, start)
    for i in start:-1:firstindex(seq)
        if f(inbounds_getindex(seq, i))
            return i
        end
    end
    return nothing
end

Base.findfirst(val, seq::BioSequence) = findnext(val, seq, 1)
Base.findlast(val, seq::BioSequence) = findprev(val, seq, lastindex(seq))


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
