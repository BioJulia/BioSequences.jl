# Kmer Iterator
# =============
#
# Iterator over all k-mers in a biological sequence.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

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
    filled::Int # This is cached values for speed
    increment::Int # This is cached values for speed
end

function each(::Type{Kmer{T,K}}, seq::Sequence, step::Integer=1) where {T,K}
    if eltype(seq) ∉ (DNA, RNA)
        throw(ArgumentError("element type must be either DNA or RNA nucleotide"))
    elseif !(0 ≤ K ≤ 32)
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

eachkmer(seq::BioSequence{A}, K::Integer, step::Integer=1) where {A<:DNAAlphabet} = each(DNAKmer{Int(K)}, seq, step)
eachkmer(seq::BioSequence{A}, K::Integer, step::Integer=1) where {A<:RNAAlphabet} = each(RNAKmer{Int(K)}, seq, step)
eachkmer(seq::ReferenceSequence, K::Integer, step::Integer=1) = each(DNAKmer{Int(K)}, seq, step)

Base.eltype(::Type{<:AbstractKmerIterator{T,S}}) where {T,S} = Tuple{Int,T}
Base.IteratorSize(::Type{<:AbstractKmerIterator{T,S}}) where {T,S} = Base.SizeUnknown()

const kmerbits = (0xff, 0x00, 0x01, 0xff,
                  0x02, 0xff, 0xff, 0xff,
                  0x03, 0xff, 0xff, 0xff,
                  0xff, 0xff, 0xff, 0xff)

function Base.iterate(it::AbstractKmerIterator{T,S}) where {T,S<:BioSequence{<:TwoBitNucs}}
    filled, i, kmer = 0, 1, UInt64(0)
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

function Base.iterate(it::EveryKmerIterator{T,S}, state::Tuple{Int,UInt64}) where {T,S<:BioSequence{<:TwoBitNucs}}
    i, kmer = state
    i += 1
    pos = i - kmersize(T) + 1
    if i > it.stop
        return nothing
    else
        nt = reinterpret(Int8, inbounds_getindex(it.seq, i))
        @inbounds val = kmerbits[nt + 1]
        kmer = kmer << 2 | val
        return (pos, T(kmer)), (i, kmer)
    end
end

function Base.iterate(it::SpacedKmerIterator{T,S}, state::Tuple{Int,UInt64}) where {T,S<:BioSequence{<:TwoBitNucs}}
    i, kmer = state
    filled = it.filled
    i += max(1, it.step - kmersize(T) + 1)
    while i ≤ it.stop
        nt = reinterpret(Int8, inbounds_getindex(it.seq, i))
        val = kmerbits[nt + 1] # inbounds
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

function Base.iterate(it::AbstractKmerIterator{T,S}) where {T,S<:BioSequence{<:FourBitNucs}}
    filled, i, kmer = 0, 1, UInt64(0)
    while i ≤ it.stop
        nt = reinterpret(Int8, inbounds_getindex(it.seq, i))
        if nt == 0xff # ambiguous
            filled = 0
        else
            filled += 1
            @inbounds val = kmerbits[nt + 1]
            kmer = kmer << 2 | val
        end
        if filled == kmersize(T)
            pos = i - kmersize(T) + 1
            return (pos, T(kmer)), (i, kmer)
        end
        i += 1
    end
    return nothing
end

function Base.iterate(it::EveryKmerIterator{T,S}, state::Tuple{Int,UInt64}) where {T,S<:BioSequence{<:FourBitNucs}}
    i, kmer = state
    i += 1
    filled = kmersize(T) - 1
    while i ≤ it.stop
        nt = reinterpret(Int8, inbounds_getindex(it.seq, i))
        if nt == 0xff # ambiguous
            filled = 0
        else
            filled += 1
            @inbounds val = kmerbits[nt + 1]
            kmer = kmer << 2 | val
        end
        if filled == kmersize(T)
            pos = i - kmersize(T) + 1
            return (pos, T(kmer)), (i, kmer)
        end
        i += 1
    end
    return nothing
end

function Base.iterate(it::SpacedKmerIterator{T,S}, state::Tuple{Int,UInt64}) where {T,S<:BioSequence{<:FourBitNucs}}
    i, kmer = state
    filled = it.filled
    i += it.increment
    while i ≤ it.stop
        nt = reinterpret(Int8, inbounds_getindex(it.seq, i))
        if nt == 0xff # ambiguous
            filled = 0
        else
            filled += 1
            @inbounds val = kmerbits[nt + 1]
            kmer = kmer << 2 | val
        end
        if filled == kmersize(T)
            pos = i - kmersize(T) + 1
            return (pos, T(kmer)), (i, kmer)
        end
        i += 1
    end
    return nothing
end
