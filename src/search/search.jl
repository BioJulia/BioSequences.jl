const DEFAULT_OVERLAP = true

struct Search{Q,I}
    query::Q
    itr::I
    overlap::Bool
end

search(query, itr; overlap = DEFAULT_OVERLAP) = Search(query, itr, overlap)

function Base.iterate(itr::Search, state=firstindex(itr.itr))
    val = findnext(itr.query, itr.itr, state)
    val === nothing && return nothing
    state = itr.overlap ? first(val) + 1 : last(val) + 1
    return val, state
end

const HasRangeEltype = Union{<:ExactSearchQuery, <:ApproximateSearchQuery, <:Regex}

Base.eltype(::Type{<:Search{Q}}) where {Q<:HasRangeEltype} = UnitRange{Int}
Base.eltype(::Type{<:Search}) = Int
Base.IteratorSize(::Type{<:Search}) = Base.SizeUnknown()

"""
    findall(pattern, sequence::BioSequence[,rng::UnitRange{Int}]; overlap::Bool=true)

Find all occurrences of `pattern` in `sequence`.

The return value is a vector of ranges of indices where the matching sequences were found.
If there are no matching sequences, the return is an empty array.

The search is restricted to the specified range when `rng` is set.

With the keyword argument `overlap` set as `true`, the start index for the next search is set to the start of the current match plus one; if set to `false`, the start index for the next search is set to the end of the current match plus one.
The default value for the keyword argument `overlap` is `true`.

See also [`ExactSearchQuery`](@ref), [`ApproximateSearchQuery`](@ref), [`PWMSearchQuery`](@ref).

# Examples
```jldoctest
julia> sequence = dna"ACACACAC"
8nt DNA Sequence:
ACACACAC

julia> findall(DNA_A, sequence)
4-element Vector{Int64}:
 1
 3
 5
 7

julia> findall(ExactSearchQuery(dna"ACAC"), sequence)
3-element Vector{UnitRange{Int64}}:
 1:4
 3:6
 5:8

julia> findall(ExactSearchQuery(dna"ACAC"), sequence; overlap=false)
2-element Vector{UnitRange{Int64}}:
 1:4
 5:8

julia> findall(ExactSearchQuery(dna"ACAC"), sequence, 2:7; overlap=false)
1-element Vector{UnitRange{Int64}}:
 3:6
```
"""
function Base.findall(pattern, sequence::BioSequence; overlap::Bool = DEFAULT_OVERLAP)
    return collect(search(pattern, sequence; overlap))
end

function Base.findall(pattern, sequence::BioSequence, rng::UnitRange{Int}; overlap::Bool = DEFAULT_OVERLAP)
    v = view(sequence, rng)
    itr = search(pattern, v; overlap)
    return map(x->parentindices(v)[1][x], itr)
end
