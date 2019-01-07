# Predicates & comparisons
# ------------------------

function Base.cmp(seq1::BioSequence, seq2::BioSequence)
    m = lastindex(seq1)
    n = lastindex(seq2)
    for i in 1:min(m, n)
        c = cmp(inbounds_getindex(seq1, i),
                inbounds_getindex(seq2, i))
        if c != 0
            return c
        end
    end
    return cmp(m, n)
end

function Base.:(==)(seq1::BioSequence, seq2::BioSequence)
    return eltype(seq1)    == eltype(seq2) &&
           length(seq1)    == length(seq2) &&
           cmp(seq1, seq2) == 0
end

Base.isless(seq1::BioSequence, seq2::BioSequence) = cmp(seq1, seq2) < 0

Base.isempty(seq::BioSequence) = length(seq) == 0

"""
    isrepetitive(seq::BioSequence, n::Integer=length(seq))

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
    ispalindromic(seq::BioSequence)

Return `true` if `seq` is a palindromic sequence; otherwise return `false`.
"""
function ispalindromic(seq::BioSequence)
    if !(eltype(seq) <: NucleicAcid)
        throw(ArgumentError("Not a nucleic acid sequence"))
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