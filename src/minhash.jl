# MinHash
# ========
#
# Functions to generate MinHash sketches of biological sequences
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

"""
MinHash Sketch type

MinHash sketches are a sorted set of kmer hashes, containing the smallest `s`
hashes for kmers of a given length. The type contains two parameters:

* .sketch: a sorted set of hashes
* .kmersize: the length of kmers used to generate the sketch
"""
struct MinHashSketch
    sketch::Vector{UInt64}
    kmersize::Int

    function MinHashSketch(sketch::Vector, kmersize::Integer)
        length(sketch) > 0 || error("Sketch cannot be empty")
        kmersize > 0 || error("Kmersize must be greater than 0")
        return new(sketch, kmersize)
    end
end

Base.getindex(s::MinHashSketch, part) = getindex(s.sketch, part)
Base.length(s::MinHashSketch) = length(s.sketch)

Base.iterate(s::MinHashSketch) = iterate(s.sketch)
Base.iterate(s::MinHashSketch, state) = iterate(s.sketch, state)

function Base.:(==)(a::MinHashSketch, b::MinHashSketch)
    return a.kmersize == b.kmersize && a.sketch == b.sketch
end


function kmerminhash!(::Type{DNAKmer{k}}, seq::MutableBioSequence, s::Integer, kmerhashes::Vector{UInt64}) where {k}
    # generate first `s` kmers
    iter = each(DNAKmer{k}, seq)
    iter_value = iterate(iter)
    while length(kmerhashes) < s && iter_value !== nothing
        (_, kmer), state = iter_value
        iter_value = iterate(iter, state)
        h = hash(canonical(kmer)) # hash lexigraphic minimum of kmer and reverse compliment of kmer
        if h ∉ kmerhashes
            push!(kmerhashes, h)
        end
    end

    sort!(kmerhashes)

    # scan `seq` to make a minhash
    while iter_value !== nothing
        (_, kmer), state = iter_value
        iter_value = iterate(iter, state)
        h = hash(canonical(kmer)) # hash lexigraphic minimum of kmer and reverse compliment of kmer
        if h < kmerhashes[end]
            i = searchsortedlast(kmerhashes, h)
            if i == 0 && h != kmerhashes[1]
                pop!(kmerhashes)
                pushfirst!(kmerhashes, h)
            elseif h != kmerhashes[i]
                pop!(kmerhashes)
                insert!(kmerhashes, i+1, h)
            end
        end
    end

    return kmerhashes
end

"""
    minhash(seq, k::Integer, s::Integer)

Generate a MinHash sketch of size `s` for kmers of length `k`.
"""
function minhash(seq::MutableBioSequence, k::Integer, s::Integer)
    kmerhashes = kmerminhash!(DNAKmer{k}, seq, s, UInt64[])
    length(kmerhashes) < s && error("failed to generate enough hashes")

    return MinHashSketch(kmerhashes, k)
end

function minhash(seqs::Vector{T}, k::Integer, s::Integer) where {T<:MutableBioSequence}
    kmerhashes = UInt64[]
    for seq in seqs
        kmerminhash!(DNAKmer{k}, seq, s, kmerhashes)
    end

    length(kmerhashes) < s && error("failed to generate enough hashes")
    return MinHashSketch(kmerhashes, k)
end

function minhash(seqs::BioCore.IO.AbstractReader, k::Integer, s::Integer)
    kmerhashes = UInt64[]
    for seq in seqs
        kmerminhash!(DNAKmer{k}, sequence(seq), s, kmerhashes)
    end

    length(kmerhashes) < s && error("failed to generate enough hashes")
    return MinHashSketch(kmerhashes, k)
end
