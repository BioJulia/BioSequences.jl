###
### Predicates & comparisons
###
### Generalised operations on biological sequence types.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

function Base.:(==)(seq1::BioSequence, seq2::BioSequence)
    seq1 === seq2 && return true
    length(seq1) == length(seq2) || return false
    @inbounds for i in eachindex(seq1)
        seq1[i] == seq2[i] || return false
    end
    return true
end

function Base.isless(seq1::BioSequence, seq2::BioSequence)
    @inbounds for i in Base.OneTo(min(length(seq1), length(seq2)))
        i1, i2 = seq1[i], seq2[i]
        isless(i1, i2) && return true
        isless(i2, i1) && return false
    end
    return isless(length(seq1), length(seq2))
end

"""
    isrepetitive(seq::BioSequence, n::Integer = length(seq))

Return `true` if and only if `seq` contains a repetitive subsequence of length `≥ n`.
"""
function isrepetitive(seq::BioSequence, n::Integer = length(seq))
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
    for i in 2:lastindex(seq)
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


"""
    ispalindromic(seq::NucSeq) -> Bool

Check if `seq` is palindromic. A palindromic sequence is identical to its
reverse-complement, so this should be equivalent to checking if
`seq == reverse_complement(seq)`.

# Examples
```jldoctest
julia> ispalindromic(dna"TGCA")
true

julia> ispalindromic(dna"TCCT")
false

julia> ispalindromic(rna"ACGGU")
false
```

Return `true` if `seq` is a palindromic sequence; otherwise return `false`.
"""
function ispalindromic(seq::BioSequence{<:NucleicAcidAlphabet})
    _ispalindromic(seq)
end

# For two-bit alphabets, all odd-length sequences are not palindromic.
function ispalindromic(seq::BioSequence{<:Union{DNAAlphabet{2}, RNAAlphabet{2}}})
    isodd(length(seq)) ? false : _ispalindromic(seq)
end

@inline function _ispalindromic(seq)
    L = lastindex(seq)
    for i in 1:cld(length(seq), 2)
        seq[i] == complement(seq[L - i + 1]) || return false
    end
    true
end

"""
    hasambiguity(seq::BioSequence)

Returns `true` if `seq` has an ambiguous symbol; otherwise return `false`.
"""
function hasambiguity(seq::BioSequence)
    for x in seq
        if isambiguous(x)
            return true
        end
    end
    return false
end
# 2 Bit specialization:
@inline hasambiguity(seq::BioSequence{<:NucleicAcidAlphabet{2}}) = false

"""
    iscanonical(seq::NucleotideSeq)

Returns `true` if `seq` is canonical.

For any sequence, there is a reverse complement, which is the same sequence, but
on the complimentary strand of DNA:

```
------->
ATCGATCG
CGATCGAT
<-------
```

!!! note
    Using the [`reverse_complement`](@ref) of a DNA sequence will give give this
    reverse complement.

Of the two sequences, the *canonical* of the two sequences is the lesser of the
two i.e. `canonical_seq < other_seq`.
"""
function iscanonical(seq::NucleotideSeq)
    i = 1
    j = lastindex(seq)
    @inbounds while i <= j
        f = seq[i]
        r = complement(seq[j])
        f < r && return true
        r < f && return false
        i += 1
        j -= 1
    end
    return true
end
