###
### A generic iterator over sequence positions that satisfy some function
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

struct ConditionIterator{F<:Function,S<:BioSequence}
    f::F
    seq::S
end

"""
    each(f::Function, seq::BioSequence) = ConditionIterator(f, seq)

Return an iterator over each element in `seq` that satisfies some function `f`.

This is similar to the [`Iterators.filter`](@ref) iterator, except that every
iteration yields a tuple of the position of the residue that satisfied `f`, and
the residue itself.
"""
each(f::Function, seq::BioSequence) = ConditionIterator(f, seq)

Base.eltype(::Type{ConditionIterator{<:Function,S}}) where {S<:BioSequence} = Tuple{Int,eltype(S)}
Base.IteratorSize(::ConditionIterator) = Base.SizeUnknown()

###
### Generic iteration for any sequence
###

Base.iterate(it::ConditionIterator) = iterate(it, findnext(it.f, it.seq, firstindex(it.seq)))
function Base.iterate(it::ConditionIterator, nextpos::Int)
    if isnothing(nextpos)
        return nothing
    else
        return (nextpos, inbounds_getindex(it.seq, nextpos)), find_next_position(it.f, it.seq, nextpos + 1)
    end
end
