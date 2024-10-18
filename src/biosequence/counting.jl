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

# These functions are BioSequences-specific because they take two arguments
isambiguous_or(x::T, y::T) where {T<:NucleicAcid} = isambiguous(x) | isambiguous(y)
isgap_or(x::T, y::T) where {T<:NucleicAcid} = isgap(x) | isgap(y)
iscertain_and(x::T, y::T) where {T<:NucleicAcid} = iscertain(x) & iscertain(y)

#BioSymbols.isambiguous(x::T, y::T) where {T<:NucleicAcid} = isambiguous(x) | isambiguous(y)
#BioSymbols.isgap(x::T, y::T) where {T<:NucleicAcid} = isgap(x) | isgap(y)
#BioSymbols.iscertain(x::T, y::T) where {T<:NucleicAcid} = iscertain(x) & iscertain(y)

Base.count(::typeof(isambiguous), seqa::S, seqb::S) where {S<:BioSequence{<:NucleicAcidAlphabet{2}}} = 0
Base.count(::typeof(isgap), seqa::S, seqb::S) where {S<:BioSequence{<:NucleicAcidAlphabet{2}}} = 0
Base.count(::typeof(iscertain), seqa::S, seqb::S) where {S<:BioSequence{<:NucleicAcidAlphabet{2}}} = min(length(seqa), length(seqb))

###
### Aliases for various uses of `count`.
###

"""
    gc_content(seq::BioSequence) -> Float64

Calculate GC content of `seq`, i.e. the number of symbols that is `DNA_C`, `DNA_G`,
`DNA_C` or `DNA_G` divided by the length of the sequence.

# Examples
```jldoctest
julia> gc_content(dna"AGCTA")
0.4

julia> gc_content(rna"UAGCGA")
0.5
```
"""
gc_content(seq::NucleotideSeq) = isempty(seq) ? 0.0 : count(isGC, seq) / length(seq)

"""
    n_ambiguous(a::BioSequence, [b::BioSequence]) -> Int

Count the number of positions where `a` (or `b`, if present) have ambigious symbols.
If `b` is given, and the length of `a` and `b` differ, look only at the indices
of the shorter sequence.

# Examples
```jldoctest
julia> n_ambiguous(dna"--TAC-WN-ACY")
3

julia> n_ambiguous(rna"UAYWW", rna"UAW")
1
```
"""
n_ambiguous(seq::BioSequence) = count(isambiguous, seq)
n_ambiguous(seqa::BioSequence, seqb::BioSequence) = count(isambiguous_or, seqa, seqb)

"""
    n_certain(a::BioSequence, [b::BioSequence]) -> Int

Count the number of positions where `a` (and `b`, if present) have certain (i.e. non-ambigous
and non-gap) symbols.
If `b` is given, and the length of `a` and `b` differ, look only at the indices
of the shorter sequence.

# Examples
```jldoctest
julia> n_certain(dna"--TAC-WN-ACY")
5

julia> n_certain(rna"UAYWW", rna"UAW")
2
```
"""
n_certain(seq::BioSequence) = count(iscertain, seq)
n_certain(seqa::BioSequence, seqb::BioSequence) = count(iscertain_and, seqa, seqb)

"""
    n_gaps(a::BioSequence, [b::BioSequence]) -> Int

Count the number of positions where `a` (or `b`, if present) have gaps.
If `b` is given, and the length of `a` and `b` differ, look only at the indices
of the shorter sequence.

# Examples
```jldoctest
julia> n_gaps(dna"--TAC-WN-ACY")
4

julia> n_gaps(dna"TC-AC-", dna"-CACG")
2
```
"""
n_gaps(seq::BioSequence) = count(isgap, seq)
n_gaps(seqa::BioSequence, seqb::BioSequence) = count(isgap_or, seqa, seqb)

"""
    mismatches(a::BioSequence, b::BioSequences) -> Int

Count the number of positions in where `a` and `b` differ.
If `b` is given, and the length of `a` and `b` differ, look only at the indices
of the shorter sequence.

# Examples
```jldoctest
julia> mismatches(dna"TAGCTA", dna"TACCTA")
1

julia> mismatches(dna"AACA", dna"AAG")
1
```
"""
mismatches(seqa::BioSequence, seqb::BioSequence) = count(!=, seqa, seqb)

"""
    mismatches(a::BioSequence, b::BioSequences) -> Int

Count the number of positions in where `a` and `b` are equal.
If `b` is given, and the length of `a` and `b` differ, look only at the indices
of the shorter sequence.

# Examples
```jldoctest
julia> matches(dna"TAGCTA", dna"TACCTA")
5

julia> matches(dna"AACA", dna"AAG")
2
```
"""
matches(seqa::BioSequence, seqb::BioSequence) = count(==, seqa, seqb)