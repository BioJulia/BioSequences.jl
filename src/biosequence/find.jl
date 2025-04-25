###
### Finders
###
###
### Finding positions meeting predicates in biological sequence types.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

Base.findfirst(f::Function, seq::BioSequence) = findnext(f, seq, firstindex(seq))
Base.findlast(f::Function, seq::BioSequence) = findprev(f, seq, lastindex(seq))

function Base.findnext(f::Function, seq::BioSequence, start::Integer)
    start = Int(start)::Int
    @boundscheck (start < 1 && throw(BoundsError(seq, start)))
    @inline _findnext(f, seq, start)
end

function Base.findprev(f::Function, seq::BioSequence, start::Integer)
    start = Int(start)::Int
    @boundscheck (start > lastindex(seq) && throw(BoundsError(seq, start)))
    @inline _findprev(f, seq, start)
end

function _findnext(f::Function, seq::BioSequence, start::Int)
    @inbounds for i in start:lastindex(seq)
        if f(seq[i])
            return i
        end
    end
    return nothing
end

function _findprev(f::Function, seq::BioSequence, start::Int)
    @inbounds for i in start:-1:firstindex(seq)
        if f(seq[i])
            return i
        end
    end
    return nothing
end

# For two-bit alphabet
# N.B: The isgap and isambiguous have identical definitions but are separate methods to avoid
# ambiguity errors.
_findnext(::typeof(isambiguous), ::BioSequence{<:NucleicAcidAlphabet{2}}, ::Int) = nothing
_findnext(::typeof(isgap), ::BioSequence{<:NucleicAcidAlphabet{2}}, ::Int) = nothing
_findprev(::typeof(isambiguous), ::BioSequence{<:NucleicAcidAlphabet{2}}, ::Int) = nothing
_findprev(::typeof(isgap), ::BioSequence{<:NucleicAcidAlphabet{2}}, ::Int) = nothing

function _findnext(::typeof(iscertain), seq::BioSequence{<:NucleicAcidAlphabet{2}}, from::Int)
    from ≤ lastindex(seq) ? from : nothing
end

function _findprev(::typeof(iscertain), seq::BioSequence{<:NucleicAcidAlphabet{2}}, from::Int)
    from ≥ 1 ? from : nothing
end

# Finding specific symbols
Base.findnext(x::DNA, seq::BioSequence{<:DNAAlphabet}, start::Integer) = findnext(isequal(x), seq, start)
Base.findnext(x::RNA, seq::BioSequence{<:RNAAlphabet}, start::Integer) = findnext(isequal(x), seq, start)
Base.findnext(x::AminoAcid, seq::BioSequence{AminoAcidAlphabet}, start::Integer) = findnext(isequal(x), seq, start)

Base.findprev(x::DNA, seq::BioSequence{<:DNAAlphabet}, start::Integer) = findprev(isequal(x), seq, start)
Base.findprev(x::RNA, seq::BioSequence{<:RNAAlphabet}, start::Integer) = findprev(isequal(x), seq, start)
Base.findprev(x::AminoAcid, seq::BioSequence{AminoAcidAlphabet}, start::Integer) = findprev(isequal(x), seq, start)

Base.findfirst(x::DNA, seq::BioSequence{<:DNAAlphabet}) = findfirst(isequal(x), seq)
Base.findfirst(x::RNA, seq::BioSequence{<:RNAAlphabet}) = findfirst(isequal(x), seq)
Base.findfirst(x::AminoAcid, seq::BioSequence{AminoAcidAlphabet}) = findfirst(isequal(x), seq)

Base.findlast(x::DNA, seq::BioSequence{<:DNAAlphabet}) = findlast(isequal(x), seq)
Base.findlast(x::RNA, seq::BioSequence{<:RNAAlphabet}) = findlast(isequal(x), seq)
Base.findlast(x::AminoAcid, seq::BioSequence{AminoAcidAlphabet}) = findlast(isequal(x), seq)

# Finding gaps
_findnext(::typeof(isgap), seq::BioSequence, i::Int) = _findnext(==(gap(eltype(seq))), seq, i)
_findprev(::typeof(isgap), seq::BioSequence, i::Int) = _findprev(==(gap(eltype(seq))), seq, i)