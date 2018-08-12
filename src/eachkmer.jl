# Kmer Iterator
# =============
#
# Iterator over all k-mers in a biological sequence.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

# Iterate through every k-mer in a nucleotide sequence
struct EachKmerIterator{T<:Kmer,S<:Sequence}
    seq::S
    step::Int
    start::Int
end

"""
    each(::Type{Kmer{T,k}}, seq::Sequence[, step=1])

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
function each(::Type{Kmer{T,K}}, seq::Sequence, step::Integer=1) where {T,K}
    if eltype(seq) ∉ (DNA, RNA)
        throw(ArgumentError("element type must be either DNA or RNA nucleotide"))
    elseif !(0 ≤ K ≤ 32)
        throw(ArgumentError("k-mer length must be between 0 and 32"))
    elseif step < 1
        throw(ArgumentError("step size must be positive"))
    end
    return EachKmerIterator{Kmer{T,K},typeof(seq)}(seq, step, 1)
end

eachkmer(seq::BioSequence{A}, K::Integer, step::Integer=1) where {A<:DNAAlphabet} = each(DNAKmer{Int(K)}, seq, step)
eachkmer(seq::BioSequence{A}, K::Integer, step::Integer=1) where {A<:RNAAlphabet} = each(RNAKmer{Int(K)}, seq, step)
eachkmer(seq::ReferenceSequence, K::Integer, step::Integer=1) = each(DNAKmer{Int(K)}, seq, step)

Base.eltype(::Type{EachKmerIterator{T,S}}) where {T,S} = Tuple{Int,T}

function Base.IteratorSize(::Type{EachKmerIterator{T,S}}) where {T,S}
    return Base.SizeUnknown()
end

Base.iterate(it::EachKmerIterator) = iterate(it, (it.start, UInt64(0)))

function Base.iterate(
        it::EachKmerIterator{T},
        state::Tuple{Int, UInt64}) where {T}

    k = kmersize(T)
    pos, kmer = state
    isok = true

    # faster path: return the next overlapping kmer if possible
    if it.step < k && pos + k - 1 ≤ lastindex(it.seq) && pos != it.start
        offset = k - it.step
        if it.step == 1
            nt = inbounds_getindex(it.seq, pos+offset)
            kmer = kmer << 2 | trailing_zeros(nt)
            isok &= iscertain(nt)
        else
            for i in 1:it.step
                nt = inbounds_getindex(it.seq, pos+i-1+offset)
                kmer = kmer << 2 | trailing_zeros(nt)
                isok &= iscertain(nt)
            end
        end
        if isok
            return (pos, T(kmer)), (pos+it.step, kmer)
        end
    end

    while pos + k - 1 ≤ lastindex(it.seq)
        kmer, ok = extract_kmer_impl(it.seq, pos, k)
        if ok
            return (pos, T(kmer)), (pos+it.step, kmer)
        end
        pos += it.step
    end

    return nothing
end

function extract_kmer_impl(seq, from, k)
    kmer::UInt64 = 0
    isok = true
    for i in 1:k
        nt = inbounds_getindex(seq, from+i-1)
        kmer = kmer << 2 | trailing_zeros(nt)
        isok &= iscertain(nt)
    end
    return kmer, isok
end
