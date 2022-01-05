
function Base.copy!(dst::BioSequence, src::BioSequence)
    resize!(dst, length(src))
    copyto!(dst, src)
end

"""
    copyto!(dst::LongSequence, src::BioSequence)

Equivalent to `copyto!(dst, 1, src, 1, length(src))`
"""
function Base.copyto!(dst::BioSequence, src::BioSequence)
    copyto!(dst, 1, src, 1, length(src))
end

"""
    copyto!(dst::BioSequence, soff, src::BioSequence, doff, N)

In-place copy `N` elements from `src` starting at `soff` to `dst`, starting at `doff`.
The length of `dst` must be greater than or equal to `N + doff - 1`.
The first N elements of `dst` are overwritten,
the other elements are left untouched. The alphabets of `src` and `dst` must be compatible.

# Examples
```
julia> seq = copyto!(dna"AACGTM", 1, dna"TAG", 1, 3)
6nt DNA Sequence:
TAGGTM

julia> copyto!(seq, 2, rna"UUUU", 1, 4)
6nt DNA Sequence:
TTTTTM
```
"""
function Base.copyto!(dst::BioSequence{A}, doff::Integer,
    src::BioSequence, soff::Integer,
    N::Integer) where {A <: Alphabet}

    @boundscheck checkbounds(dst, doff:doff+N-1)
    @boundscheck checkbounds(src, soff:soff+N-1)

    for i in 0:N-1
        dst[doff + i] = src[soff + i]
    end
    
    return dst
end