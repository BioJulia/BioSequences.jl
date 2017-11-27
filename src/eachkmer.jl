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

abstract type AbstractKmerIterator{T,S} end

struct EveryKmerIterator{T<:Kmer,S<:Sequence} <: AbstractKmerIterator{T,S}
    seq::S
    start::Int
    stop::Int
end

struct SpacedKmerIterator{T<:Kmer,S<:Sequence} <: AbstractKmerIterator{T,S}
    seq::S
    start::Int
    step::Int
    stop::Int
    filled::Int # This is cached for speed
    increment::Int # This is cached for speed
end

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
function each(::Type{Kmer{T,K}}, seq::BioSequence, step::Integer=1) where {T,K}
    if eltype(seq) ∉ (DNA, RNA)
        throw(ArgumentError("element type must be either DNA or RNA nucleotide"))
    elseif !(1 ≤ K ≤ 32)
        throw(ArgumentError("k-mer length must be between 0 and 32"))
    elseif step < 1
        throw(ArgumentError("step size must be positive"))
    end
    if step == 1
        return EveryKmerIterator{Kmer{T,K},typeof(seq)}(seq, 1, lastindex(seq))
    else
        filled = max(0, K - step)
        increment = max(1, step - K + 1)
        return SpacedKmerIterator{Kmer{T,K},typeof(seq)}(seq, 1, step, lastindex(seq), filled, increment)
    end
end

function eachkmer(seq::MutableBioSequence{A}, K::Integer, step::Integer=1) where {A<:DNAAlphabet}
    Base.depwarn("eachkmer is depreceated: type instability means it is too slow. Please use each(::Type{Kmer{T,K}}, seq, step) instead", :eachkmer)
    return each(DNAKmer{Int(K)}, seq, step)
end
function eachkmer(seq::MutableBioSequence{A}, K::Integer, step::Integer=1) where {A<:RNAAlphabet}
    Base.depwarn("eachkmer is depreceated: type instability means it is too slow. Please use each(::Type{Kmer{T,K}}, seq, step) instead", :eachkmer)
    return each(RNAKmer{Int(K)}, seq, step)
end
function eachkmer(seq::ReferenceSequence, K::Integer, step::Integer=1)
    Base.depwarn("eachkmer is depreceated: type instability means it is too slow. Please use each(::Type{Kmer{T,K}}, seq, step) instead", :eachkmer)
    return each(DNAKmer{Int(K)}, seq, step)
end

Base.eltype(::Type{<:AbstractKmerIterator{T,S}}) where {T,S} = Tuple{Int,T}
Base.IteratorSize(::Type{<:AbstractKmerIterator{T,S}}
) where {T,S<:Union{ReferenceSequence, BioSequence{<:FourBitNucs}}} = Base.SizeUnknown()
Base.IteratorSize(::Type{<:AbstractKmerIterator{T,S}}
) where {T,S<:BioSequence{<:TwoBitNucs}} = Base.HasLength()

function Base.length(it::AbstractKmerIterator{T,S}) where {T,S<:BioSequence{<:TwoBitNucs}}
    return max(0, fld(it.stop - it.start + 1 - kmersize(T), step(it)) + 1)
end

Base.step(x::EveryKmerIterator) = 1
Base.step(x::SpacedKmerIterator) = x.step

# A nucleotide with bitvalue B has kmer-bitvalue kmerbits[B+1].
# ambiguous nucleotides have no kmervalue, here set to 0xff
const kmerbits = (0xff, 0x00, 0x01, 0xff,
                  0x02, 0xff, 0xff, 0xff,
                  0x03, 0xff, 0xff, 0xff,
                  0xff, 0xff, 0xff, 0xff)

# Initializer for TwoBitNucs
function Base.iterate(it::AbstractKmerIterator{T,S}) where {T,S<:BioSequence{<:TwoBitNucs}}
    filled, i, kmer = 0, it.start, UInt64(0)

    while i ≤ it.stop
        nt = reinterpret(Int8, inbounds_getindex(it.seq, i))
        @inbounds val = kmerbits[nt + 1]
        kmer = kmer << 2 | val
        filled += 1
        if filled == kmersize(T)
            return (1, T(kmer)), (i, kmer)
        end
        i += 1
    end
    return nothing
end

function Base.iterate(it::EveryKmerIterator{T,S}, state
    ) where {T,S<:BioSequence{<:TwoBitNucs}}
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

function Base.iterate(it::SpacedKmerIterator{T,S}, state
    ) where {T,S<:BioSequence{<:TwoBitNucs}}
    i, kmer = state
    filled = it.filled
    i += it.increment

    while i ≤ it.stop
        nt = reinterpret(Int8, inbounds_getindex(it.seq, i))
        @inbounds val = kmerbits[nt + 1]
        kmer = kmer << 2 | val
        filled += 1
        if filled == kmersize(T)
            pos = i - kmersize(T) + 1
            return (pos, T(kmer)), (i, kmer)
        end
        i += 1
    end
    return nothing
end

function Base.iterate(it::EveryKmerIterator{T,S}, state=(it.start-1,1,UInt64(0))
    ) where {T,S<:Union{ReferenceSequence, BioSequence{<:FourBitNucs}}}
    i, filled, kmer = state
    i += 1
    filled -= 1

    while i ≤ it.stop
        nt = reinterpret(Int8, inbounds_getindex(it.seq, i))
        @inbounds val = kmerbits[nt + 1]
        kmer = kmer << 2 | val
        filled = ifelse(val == 0xff, 0, filled+1)

        if filled == kmersize(T)
            pos = i - kmersize(T) + 1
            return (pos, T(kmer)), (i, kmersize(T), kmer)
        end
        i += 1
    end
    return nothing
end

@inline function Base.iterate(it::SpacedKmerIterator{T,S}, state=(it.start-it.increment, 1, 0, UInt64(0))
    ) where {T,S<:Union{ReferenceSequence, BioSequence{<:FourBitNucs}}}
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
        if filled == kmersize(T)
            state = (i, i - kmersize(T) + 1 + it.step, it.filled, kmer)
            return (pos, T(kmer)), state
        end
        i += 1
    end
    return nothing
end
