# Transformations 
# ===============
#
# Methods that manipulate and change a biological sequence.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md



"""
    empty!(seq)

Completely empty a biological sequence `seq` of nucleotides.
"""
Base.empty!(seq::BioSequence) = resize!(seq, 0)

"""
    push!(seq, x)

Append a biological symbol `x` to a biological sequence `seq`.
"""
function Base.push!(seq::BioSequence, x)
    if eltype(seq) !== typeof(x)
        throw(ArgumentError(string("Can't push a ", typeof(x), " onto a ", eltype(seq), " sequence.")))
    end
    resize!(seq, length(seq) + 1)
    unsafe_setindex!(seq, x, lastindex(seq))
    return seq
end

"""
    pop!(seq::BioSequence)

Remove the symbol from the end of a biological sequence `seq` and return it.
Returns a variable of `eltype(seq)`.
"""
function Base.pop!(seq::BioSequence)
    if isempty(seq)
        throw(ArgumentError("sequence must be non-empty"))
    end
    @inbounds x = seq[end]
    deleteat!(seq, lastindex(seq))
    return x
end

"""
    insert!(seq, i, x)

Insert a biological symbol `x` into a biological sequence `seq`, at the given
index `i`.
"""
function Base.insert!(seq::BioSequence, i::Integer, x)
    if eltype(seq) !== typeof(x)
        throw(ArgumentError(string("Can't insert a ", typeof(x), " into a ", eltype(seq), " sequence.")))
    end
    checkbounds(seq, i)
    resize!(seq, length(seq) + 1)
    copyto!(seq, i + 1, seq, i, lastindex(seq) - i)
    unsafe_setindex!(seq, x, i)
    return seq
end