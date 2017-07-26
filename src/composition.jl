# Composition
# ===========
#
# Sequence composition counter.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

struct Composition{T} <: Associative{T,Int}
    counts::Vector{Int}
end

function Composition(seq::BioSequence{A}) where {A<:Union{DNAAlphabet,RNAAlphabet}}
    counts = zeros(Int, 16)
    @inbounds for x in seq
        counts[reinterpret(UInt8, x) + 1] += 1
    end
    return Composition{eltype(A)}(counts)
end

function Composition(seq::ReferenceSequence)
    counts = zeros(Int, 16)
    @inbounds for x in seq
        counts[reinterpret(UInt8, x) + 1] += 1
    end
    return Composition{DNA}(counts)
end

function Composition(kmer::Kmer{T,k}) where {T,k}
    counts = zeros(Int, 16)
    @inbounds begin
        counts[Int(DNA_A)+1] = count_a(kmer)
        counts[Int(DNA_C)+1] = count_c(kmer)
        counts[Int(DNA_G)+1] = count_g(kmer)
        counts[Int(DNA_T)+1] = count_t(kmer)  # U when T == RNA
    end
    return Composition{T}(counts)
end

function Composition(seq::AminoAcidSequence)
    counts = zeros(Int, length(alphabet(AminoAcid)))
    @inbounds for x in seq
        counts[reinterpret(UInt8, x) + 1] += 1
    end
    return Composition{AminoAcid}(counts)
end

function Composition(iter::EachKmerIterator{T}) where {T<:Kmer}
    counts = zeros(Int, 4^kmersize(T))
    for (_, x) in iter
        counts[convert(UInt64, x) + 1] += 1
    end
    return Composition{T}(counts)
end

"""
    composition(seq | kmer_iter)

Calculate composition of biological symbols in `seq` or k-mers in `kmer_iter`.
"""
function composition(iter::Union{Sequence,EachKmerIterator})
    return Composition(iter)
end

function Base.:(==){T}(x::Composition{T}, y::Composition{T})
    return x.counts == y.counts
end

function Base.length(comp::Composition)
    return length(comp.counts)
end

function Base.start(comp::Composition)
    return UInt64(0)
end

function Base.done(comp::Composition, i)
    return i ≥ length(comp.counts)
end

function Base.next(comp::Composition, i)
    key = convert(keytype(comp), i)
    count = comp.counts[i + 1]
    return (key => count), i + 1
end

function Base.getindex(comp::Composition{T}, x) where {T}
    i = convert(UInt64, convert(T, x))
    if !(0 ≤ i < endof(comp.counts))
        throw(KeyError(x))
    end
    return comp.counts[i+1]
end

function Base.copy(comp::Composition{T}) where {T}
    return Composition{T}(copy(comp.counts))
end

function Base.merge(comp::Composition{T}, other::Composition{T}) where {T}
    return merge!(copy(comp), other)
end

function Base.merge!(comp::Composition{T}, other::Composition{T}) where {T}
    @assert length(comp.counts) == length(other.counts)
    for i in 1:endof(comp.counts)
        comp.counts[i] += other.counts[i]
    end
    return comp
end

function Base.summary(::Composition{T}) where {T}
    if T == DNA
        return "DNA Composition"
    elseif T == RNA
        return "RNA Composition"
    elseif T == AminoAcid
        return "Amino Acid Composition"
    else
        return string(T, " Composition")
    end
end
