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
        throw(ArgumentError("Can't push a $(typeof(x)) onto a sequence of $(eltype(seq))"))
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