###
### Finders
###
###
### Finding positions meeting predicates in biological sequence types.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

function Base.findnext(f::Function, seq::BioSequence, start::Integer)
    start > lastindex(seq) && return nothing
    checkbounds(seq, start)
    @inbounds for i in start:lastindex(seq)
        if f(seq[i])
            return i
        end
    end
    return nothing
end

# No ambiguous sites can exist in a nucleic acid sequence using the two-bit alphabet.
Base.findnext(::typeof(isambiguous), seq::BioSequence{<:NucleicAcidAlphabet{2}}, from::Integer) = nothing

function Base.findprev(f::Function, seq::BioSequence, start::Integer)
    start < firstindex(seq) && return nothing
    checkbounds(seq, start)
    for i in start:-1:firstindex(seq)
        if f(seq[i])
            return i
        end
    end
    return nothing
end

Base.findfirst(f::Function, seq::BioSequence) = findnext(f, seq, 1)
Base.findlast(f::Function, seq::BioSequence) = findprev(f, seq, lastindex(seq))
