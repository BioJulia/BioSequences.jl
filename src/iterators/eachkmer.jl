# Kmer Iterator
# =============
#
# Iterator over all k-mers in a biological sequence.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

# Note about the variable names:

# it.step is the distance from start of one kmer to start of next
#
# filled is the number of nucleotides in a kmer that has the correct value set.
# e.g. when moving from xxxxxxx to yyyyyyy, filled goes from K to 3.
# when moving from xxxxxxx to zzzzzzz, filled goes from K to 0.
#
# increment is how far the index jumps ahead when going to the next kmer.
# for close kmers where no jump is possible, it's 1, else it can be arbitrary high.
#
# For a kmeriterator that first emits xxxxxxx, then zzzzzzz:
#
#       |--------- step ---------|
#       xxxxxxx                  zzzzzzz
# -------------------------------------------------------------
#           yyyyyyy
#             |---- increment ---|
#
# The state returned at each iteration is the state upon return, not the state
# needed for the following iteration.

abstract type AbstractMerIterator{T,S,C} end

struct EveryMerIterator{T<:AbstractMer,S<:BioSequence,C} <: AbstractMerIterator{T,S,C}
    seq::S
    start::Int
    stop::Int
end

struct SpacedMerIterator{T<:AbstractMer,S<:BioSequence,C} <: AbstractMerIterator{T,S,C}
    seq::S
    start::Int
    step::Int
    stop::Int
    filled::Int # This is cached for speed
    increment::Int # This is cached for speed
end

#=
struct Canonical{I<:AbstractKmerIterator}
    it::I
end
@inline Base.length(x::Canonical) = length(x.it)
@inline Base.eltype(x::Canonical) = eltype(x.it)
@inline Base.IteratorSize(x::Canonical) = Base.IteratorSize(x.it)

@inline function Base.iterate(x::Canonical)
    i = iterate(x.it)
    i === nothing && return nothing
    @inbounds pos, mer = i[1]
    @inbounds return (pos, canonical(mer)), i[2]
end

@inline function Base.iterate(x::Canonical, state)
    i = iterate(x.it, state)
    i === nothing && return nothing
    @inbounds pos, mer = i[1]
    @inbounds return (pos, canonical(mer)), i[2]
end
=#

"""
    each(::Type{Kmer{T,k}}, seq::BioSequence[, step=1])

Initialize an iterator over all k-mers in a sequence `seq` skipping ambiguous
nucleotides without changing the reading frame.

# Arguments
* `Kmer{T,k}`: k-mer type to enumerate.
* `seq`: a nucleotide sequence.
* `step=1`: the number of positions between iterated k-mers

# Examples
```
# iterate over DNA codons
for (pos, codon) in each(DNAKmer{3}, dna"ATCCTANAGNTACT", 3)
    @show pos, codon
end
```
"""
#function each(::Type{Skipmer{U, A, M, M, K}}, seq::BioSequence, step::Integer = 1) where {U, A, M, K}
function each(::Type{T}, seq::BioSequence{A}, step::Integer = 1) where {A<:NucleicAcidAlphabet{2},K,T<:AbstractMer{A,K}}
    if !(1 ≤ K ≤ capacity(T))
        throw(ArgumentError("k-mer length must be between 0 and $(capacity(T))"))
    elseif step < 1
        throw(ArgumentError("step size must be positive"))
    end
    if step == 1
        return EveryMerIterator{T,typeof(seq),false}(seq, 1, lastindex(seq))
    else
        filled = max(0, K - step)
        increment = max(1, step - K + 1)
        return SpacedMerIterator{T,typeof(seq), false}(seq, 1, step, lastindex(seq), filled, increment)
    end
end

function eachcanonical(::Type{T}, seq::BioSequence{A}, step::Integer = 1) where {A<:NucleicAcidAlphabet{2},K,T<:AbstractMer{A,K}}
    if !(1 ≤ K ≤ capacity(T))
        throw(ArgumentError("k-mer length must be between 0 and $(capacity(T))"))
    elseif step < 1
        throw(ArgumentError("step size must be positive"))
    end
    if step == 1
        return EveryMerIterator{T,typeof(seq),true}(seq, 1, lastindex(seq))
    else
        filled = max(0, K - step)
        increment = max(1, step - K + 1)
        return SpacedMerIterator{T,typeof(seq),true}(seq, 1, step, lastindex(seq), filled, increment)
    end
end

@inline Base.eltype(::Type{<:AbstractMerIterator{T,S,C}}) where {T,S,C} = Tuple{Int,T}
@inline Base.IteratorSize(::Type{<:AbstractMerIterator{T,S,C}}
) where {T,S<:Union{ReferenceSequence, LongSequence{<:NucleicAcidAlphabet{4}}},C} = Base.SizeUnknown()
@inline Base.IteratorSize(::Type{<:AbstractMerIterator{T,S,C}}
) where {T,S<:LongSequence{<:NucleicAcidAlphabet{2}},C} = Base.HasLength()

@inline function Base.length(it::AbstractMerIterator{T,S,C}) where {T,S<:LongSequence{<:NucleicAcidAlphabet{2}},C}
    return max(0, fld(it.stop - it.start + 1 - K(T), step(it)) + 1)
end

Base.step(x::EveryMerIterator) = 1
Base.step(x::SpacedMerIterator) = x.step

# A nucleotide with bitvalue B has kmer-bitvalue kmerbits[B+1].
# ambiguous nucleotides have no kmervalue, here set to 0xff
const kmerbits = (0xff, 0x00, 0x01, 0xff,
                  0x02, 0xff, 0xff, 0xff,
                  0x03, 0xff, 0xff, 0xff,
                  0xff, 0xff, 0xff, 0xff)

# Initializers for two-bit nucleic acid alphabets...
@inline function Base.iterate(it::AbstractMerIterator{T,S,false}) where {T,S<:LongSequence{<:NucleicAcidAlphabet{2}}}
    filled, i, kmer = 0, it.start, zero(encoded_data_type(T))

    while i ≤ it.stop
        nt = reinterpret(Int8, inbounds_getindex(it.seq, i))
        @inbounds val = kmerbits[nt + 1]
        kmer = kmer << 2 | val
        filled += 1
        if filled == K(T)
            return (1, T(kmer)), (i, kmer)
        end
        i += 1
    end
    return nothing
end

@inline function Base.iterate(it::AbstractMerIterator{T,S,true}) where {T,S<:LongSequence{<:NucleicAcidAlphabet{2}}}
    UT = encoded_data_type(T)
    filled, i = 0, it.start
    fwkmer = rvkmer = zero(UT)

    while i ≤ it.stop
        nt = reinterpret(UInt8, inbounds_getindex(it.seq, i))
        @inbounds fbits = UT(kmerbits[nt + 1])
        rbits = ~fbits & typeof(fbits)(0x03)
        fwkmer = (fwkmer << 0x02 | fbits)
        rvkmer = (rvkmer >> 0x02) | (UT(rbits) << unsigned(offset(T, 1)))
        filled += 1
        if filled == K(T)
            return (1, min(T(fwkmer), T(rvkmer))), (i, fwkmer, rvkmer)
        end
        i += 1
    end
    return nothing
end

@inline function Base.iterate(it::EveryMerIterator{T,S,false}, state
    ) where {T,S<:LongSequence{<:NucleicAcidAlphabet{2}}}
    i, kmer = state
    i += 1

    if i > it.stop
        return nothing
    else
        nt = reinterpret(Int8, inbounds_getindex(it.seq, i))
        @inbounds val = kmerbits[nt + 1]
        kmer = kmer << 2 | val
        pos = i - kmersize(T) + 1
        return (pos, T(kmer)), (i, kmer)
    end
end

@inline function Base.iterate(it::EveryMerIterator{T,S,true}, state) where {T, S<:LongSequence{<:NucleicAcidAlphabet{2}}}
    UT = encoded_data_eltype(T)
    i, fwkmer, rvkmer = state
    i += 1
    
    if i > it.stop
        return nothing
    else
        nt = reinterpret(UInt8, inbounds_getindex(it.seq, i))
        @inbounds fbits = UT(kmerbits[nt + 1])
        rbits = ~fbits & typeof(fbits)(0x03)
        fwkmer = (fwkmer << 0x02 | fbits)
        rvkmer = (rvkmer >> 0x02) | (UT(rbits) << unsigned(offset(T, 1)))
        pos = i - kmersize(T) + 1
        
        return (pos, min(T(fwkmer), T(rvkmer))), (i, fwkmer, rvkmer)
    end    
end

@inline function Base.iterate(it::SpacedMerIterator{T,S,false}, state
    ) where {T,S<:LongSequence{<:NucleicAcidAlphabet{2}}}
    i, kmer = state
    filled = it.filled
    i += it.increment

    while i ≤ it.stop
        nt = reinterpret(Int8, inbounds_getindex(it.seq, i))
        @inbounds val = kmerbits[nt + 1]
        kmer = kmer << 2 | val
        filled += 1
        if filled == K(T)
            pos = i - K(T) + 1
            return (pos, T(kmer)), (i, kmer)
        end
        i += 1
    end
    return nothing
end

@inline function Base.iterate(it::SpacedMerIterator{T,S,true}, state
    ) where {T,S<:LongSequence{<:NucleicAcidAlphabet{2}}}
    UT = encoded_data_type(T)
    i, fkmer, rkmer = state
    filled = it.filled
    i += it.increment

    while i ≤ it.stop
        nt = reinterpret(Int8, inbounds_getindex(it.seq, i))
        @inbounds val = kmerbits[nt + 1]
        rbits = ~val & typeof(val)(0x03)
        fkmer = fkmer << 2 | val
        rkmer = rkmer >> 2 | (UT(rbits) << unsigned(offset(T, 1)))
        filled += 1
        if filled == K(T)
            pos = i - K(T) + 1
            return (pos, min(T(fkmer), T(rkmer))), (i, fkmer, rkmer)
        end
        i += 1
    end
    return nothing
end

@inline function Base.iterate(it::EveryMerIterator{T,S,false}, state=(it.start-1,1,UInt64(0))
    ) where {T,S<:Union{ReferenceSequence, LongSequence{<:NucleicAcidAlphabet{4}}}}
    i, filled, kmer = state
    i += 1
    filled -= 1

    while i ≤ it.stop
        nt = reinterpret(Int8, inbounds_getindex(it.seq, i))
        @inbounds val = kmerbits[nt + 1]
        kmer = kmer << 2 | val
        filled = ifelse(val == 0xff, 0, filled + 1)

        if filled == K(T)
            pos = i - K(T) + 1
            return (pos, T(kmer)), (i, K(T), kmer)
        end
        i += 1
    end
    return nothing
end

@inline function Base.iterate(it::EveryMerIterator{T,S,true}, state=(it.start-1,1,encoded_data_type(T)(0),encoded_data_type(T)(0))
    ) where {T,S<:Union{ReferenceSequence, LongSequence{<:NucleicAcidAlphabet{4}}}}
    
    UT = encoded_data_type(T)
    i, filled, fkmer, rkmer = state
    i += 1
    filled -= 1

    while i ≤ it.stop
        nt = reinterpret(Int8, inbounds_getindex(it.seq, i))
        @inbounds fbits = UT(kmerbits[nt + 1])
        rbits = ~fbits & typeof(fbits)(0x03)
        fkmer = fkmer << 0x02 | fbits
        rkmer = (rkmer >> 0x02) | (UT(rbits) << unsigned(offset(T, 1)))
        filled = ifelse(fbits == 0xff, 0, filled + 1)
        if filled == K(T)
            pos = i - K(T) + 1
            return (pos, min(T(fkmer), T(rkmer))), (i, K(T), fkmer, rkmer)
        end
        i += 1
    end
    return nothing
end

@inline function Base.iterate(it::SpacedMerIterator{T,S,false}, state=(it.start-it.increment, 1, 0, zero(encoded_data_type(T)))
    ) where {T,S<:Union{ReferenceSequence, LongSequence{<:NucleicAcidAlphabet{4}}}}
    i, pos, filled, kmer = state
    i += it.increment

    while i ≤ it.stop
        nt = reinterpret(Int8, inbounds_getindex(it.seq, i))
        @inbounds val = kmerbits[nt + 1]
        if val == 0xff # ambiguous
            filled = 0
            # Find the beginning of next possible kmer after i
            pos = i + it.step - Core.Intrinsics.urem_int(i-pos, it.step)
            i = pos - 1
        else
            filled += 1
            kmer = kmer << 2 | val
        end
        if filled == K(T)
            state = (i, i - K(T) + 1 + it.step, it.filled, kmer)
            return (pos, T(kmer)), state
        end
        i += 1
    end
    return nothing
end

@inline function Base.iterate(it::SpacedMerIterator{T,S,true}, state=(it.start-it.increment, 1, 0, zero(encoded_data_type(T)), zero(encoded_data_type(T)))
    ) where {T,S<:Union{ReferenceSequence, LongSequence{<:NucleicAcidAlphabet{4}}}}
    UT = encoded_data_eltype(T)
    i, pos, filled, fwkmer, rvkmer = state
    i += it.increment

    while i ≤ it.stop
        nt = reinterpret(Int8, inbounds_getindex(it.seq, i))
        @inbounds val = kmerbits[nt + 1]
        rbits = ~val & typeof(val)(0x03)
        if val == 0xff # ambiguous
            filled = 0
            # Find the beginning of next possible kmer after i
            pos = i + it.step - Core.Intrinsics.urem_int(i-pos, it.step)
            i = pos - 1
        else
            filled += 1
            fwkmer = fwkmer << 2 | val
            rvkmer = (rvkmer >> 0x02) | UT(rbits) << unsigned(offset(T, 1))
        end
        if filled == K(T)
            state = (i, i - K(T) + 1 + it.step, it.filled, fwkmer, rvkmer)
            return (pos, min(T(fwkmer), T(rvkmer))), state
        end
        i += 1
    end
    return nothing
end
