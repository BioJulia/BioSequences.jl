###
### Mer Iteration
###
###
### Iterator over all k-mers in a biological sequence.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

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

struct MerIterResult{T<:AbstractMer}
    position::Int
    fw::T
    bw::T
end

@inline function Base.iterate(x::MerIterResult, state = 1)
    if 1 ≤ state ≤ 3
        return getfield(x, state), state + 1
    else
        return nothing
    end
end

@inline position(x::MerIterResult) = x.position
@inline fwmer(x::MerIterResult) = x.fw
@inline bwmer(x::MerIterResult) = x.bw
@inline canonical(x::MerIterResult) = min(fwmer(x), bwmer(x))

function Base.show(io::IO, ::MIME"text/plain", x::MerIterResult)
    println(io, "Mer iteration result:")
    println(io, "Position: ", x.position)
    print(io, "Forward: ")
    showcompact(io, fwmer(x))
    print(io, "\nBackward: ")
    showcompact(io, bwmer(x))
    print(io, '\n')
end

abstract type AbstractMerIterator{T,S} end

struct EveryMerIterator{T<:AbstractMer,S<:BioSequence} <: AbstractMerIterator{T,S}
    seq::S
    start::Int
    stop::Int
end

struct SpacedMerIterator{T<:AbstractMer,S<:BioSequence} <: AbstractMerIterator{T,S}
    seq::S
    start::Int
    step::Int
    stop::Int
    filled::Int # This is cached for speed
    increment::Int # This is cached for speed
end

"""
    each(::Type{T}, seq::BioSequence) where {T<:AbstractMer}

Initialize an iterator over all overlapping k-mers in a sequence `seq` skipping
ambiguous nucleotides without changing the reading frame.
"""
function each(::Type{T}, seq::BioSequence) where {T<:AbstractMer}
    if eltype(seq) ∉ (DNA, RNA)
        throw(ArgumentError("element type must be either DNA or RNA nucleotide"))
    elseif !(1 ≤ ksize(T) ≤ capacity(T))
        throw(ArgumentError("k-mer length must be between 0 and $(capacity(T))"))
    end
    return EveryMerIterator{T,typeof(seq)}(seq, 1, lastindex(seq))
end

"""
    each(::Type{T}, seq::BioSequence, step::Integer) where {T<:AbstractMer}

Initialize an iterator over k-mers separated by a `step` parameter, in a
sequence `seq` skipping ambiguous nucleotides without changing the reading frame.
"""
function each(::Type{T}, seq::BioSequence, step::Integer) where {T<:AbstractMer}
    if eltype(seq) ∉ (DNA, RNA)
        throw(ArgumentError("element type must be either DNA or RNA nucleotide"))
    elseif !(1 ≤ ksize(T) ≤ capacity(T))
        throw(ArgumentError("k-mer length must be between 0 and $(capacity(T))"))
    elseif step < 1
        throw(ArgumentError("step size must be positive"))
    end
    if step == 1
        return EveryMerIterator{T,typeof(seq)}(seq, 1, lastindex(seq))
    else
        filled = max(0, ksize(T) - step)
        increment = max(1, step - ksize(T) + 1)
        return SpacedMerIterator{T,typeof(seq)}(seq, 1, step, lastindex(seq), filled, increment)
    end
end

@inline Base.eltype(::Type{<:AbstractMerIterator{T,S}}) where {T,S} = MerIterResult{T}
@inline Base.IteratorSize(::Type{<:AbstractMerIterator{T,S}}) where {T,S<:LongSequence{<:NucleicAcidAlphabet{4}}} = Base.SizeUnknown()
@inline Base.IteratorSize(::Type{<:AbstractMerIterator{T,S}}) where {T,S<:LongSequence{<:NucleicAcidAlphabet{2}}} = Base.HasLength()

@inline function Base.length(it::AbstractMerIterator{T,S}) where {T,S<:LongSequence{<:NucleicAcidAlphabet{2}}}
    return max(0, fld(it.stop - it.start + 1 - ksize(T), step(it)) + 1)
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
@inline function Base.iterate(it::AbstractMerIterator{T,S}) where {T,S<:LongSequence{<:NucleicAcidAlphabet{2}}}
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
        if filled == ksize(T)
            return MerIterResult(1, T(fwkmer), T(rvkmer)), (i, fwkmer, rvkmer)
        end
        i += 1
    end
    return nothing
end

@inline function Base.iterate(it::EveryMerIterator{T,S}, state) where {T,S<:LongSequence{<:NucleicAcidAlphabet{2}}}
    UT = encoded_data_type(T)
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
        pos = i - ksize(T) + 1
        
        return MerIterResult(pos, T(fwkmer), T(rvkmer)), (i, fwkmer, rvkmer)
    end    
end

@inline function Base.iterate(it::SpacedMerIterator{T,S}, state
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
        if filled == ksize(T)
            pos = i - ksize(T) + 1
            return MerIterResult(pos, T(fkmer), T(rkmer)), (i, fkmer, rkmer)
        end
        i += 1
    end
    return nothing
end

@inline function Base.iterate(it::EveryMerIterator{T,S}, state=(it.start-1,1,encoded_data_type(T)(0),encoded_data_type(T)(0))
    ) where {T,S<:LongSequence{<:NucleicAcidAlphabet{4}}}
    
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
        if filled == ksize(T)
            pos = i - ksize(T) + 1
            return MerIterResult(pos, T(fkmer), T(rkmer)), (i, ksize(T), fkmer, rkmer)
        end
        i += 1
    end
    return nothing
end

@inline function Base.iterate(it::SpacedMerIterator{T,S}, state=(it.start-it.increment, 1, 0, zero(encoded_data_type(T)), zero(encoded_data_type(T)))
    ) where {T,S<:LongSequence{<:NucleicAcidAlphabet{4}}}
    UT = encoded_data_type(T)
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
        if filled == ksize(T)
            state = (i, i - ksize(T) + 1 + it.step, it.filled, fwkmer, rvkmer)
            return MerIterResult(pos, T(fwkmer), T(rvkmer)), state
        end
        i += 1
    end
    return nothing
end
