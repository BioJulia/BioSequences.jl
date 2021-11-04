# Exact Search
# ============
#
# Exact sequence search tools.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

# Reference Implementation
# ------------------------
#
# `searchindex` and `rsearchindex` should work exactly the same way as the
# following reference implementation:
#=
# Return the index of the first occurrence of `pat` in `seq[start:stop]`.
function searchindex(seq, pat, start=1, stop=lastindex(seq))
    m = length(pat)
    n = length(seq)
    for s in max(start-1, 0):min(stop, n)-m
        if occurs_with_shift(pat, seq, s)
            return s+1  # found
        end
    end
    return 0  # not found
end

# Return the index of the last occurrence of `pat` in `seq[stop:start]`.
function rsearchindex(seq, pat, start=lastindex(seq), stop=1)
    n = length(seq)
    m = length(pat)
    for s in min(start, n)-m:-1:max(stop-1, 0)
        if occurs_with_shift(pat, seq, s)
            return s+1  # found
        end
    end
    return 0  # not found
end

function occurs_with_shift(pat, seq, shift)
    for i in 1:lastindex(pat)
        if !iscompatible(pat[i], seq[shift+i])
            return false
        end
    end
    return true
end
=#


# Actual Implementation
# ---------------------


# Forward
# -------

# This algorithm is borrowed from base/strings/search.jl, which looks similar to
# Sunday's Quick Search algorithm.


# Backward
# --------

