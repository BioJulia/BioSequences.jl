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

# Thought of more clever ways of doing this but for 3 lines I think this is fine
# for exact symbol finding. Plus it will error is you do say findfirst(DNA_A, rna"AUG").
# Whilst it is concieveable to want the DNA to be converted to RNA for search,
# that would be an implicit behaviour some may no expect. Probably best to error
# and cause a fuss if one was trying to search for RNA letters in a container type
# that does not contain them!
Base.findfirst(x::DNA, seq::BioSequence{<:DNAAlphabet}) = Base.findfirst(isequal(x), seq)
Base.findfirst(x::RNA, seq::BioSequence{<:RNAAlphabet}) = Base.findfirst(isequal(x), seq)
Base.findfirst(x::AminoAcid, seq::BioSequence{AminoAcidAlphabet}) = Base.findfirst(isequal(x), seq)