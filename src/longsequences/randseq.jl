# Random Sequence Generator
# =========================
#
# Random sequence generator.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

import Random: Sampler

"""
    SamplerUniform{T}

Weighted sampler of type T. Instantiate with a collection of eltype T containing
the elements to sample.

# Examples
```
julia> sp = SamplerUniform(rna"ACGU");
```
"""
struct SamplerUniform{T} <: Sampler{T}
    elems::Vector{T}

    function SamplerUniform{T}(elems) where {T}
        elemsvector = convert(Vector{T}, collect(elems))
        return new(elemsvector)
    end
end

SamplerUniform(elems) = SamplerUniform{eltype(elems)}(elems)
Base.rand(rng::AbstractRNG, sp::SamplerUniform) = rand(sp.elems)
const DefaultAASampler = SamplerUniform(aa"ACDEFGHIKLMNPQRSTVWY")

"""
    SamplerWeighted{T}

Weighted sampler of type T. Instantiate with a collection of eltype T containing
the elements to sample, and an orderen collection of probabilities to sample
each element except the last. The last probability is the remaining probability
up to 1.

# Examples
```
julia> sp = SamplerWeighted(rna"ACGUN", fill(0.2475, 4));
```
"""
struct SamplerWeighted{T} <: Sampler{T}
    elems::Vector{T}
    probs::Vector{Float64}

    function SamplerWeighted{T}(elems, probs) where {T}
        probsum = sum(probs)
        if probsum ≥ 1.0
            throw(ArgumentError("sum of probabilties cannot exceed 1.0"))
        elseif length(elems) != length(probs) + 1
            throw(ArgumentError("length of elems must be length of probs + 1"))
        elseif minimum(probs) < 0.0
            throw(ArgumentError("probabilities must be non-negative"))
        end
        # Even with float weirdness,  we can guarantee x + (1.0 - x) == 1.0,
        # when 0 ≤ x ≤ 1, as there's no exponent, and it works like int addition
        elemsvector = convert(Vector{T}, collect(elems))
        probsvector = push!(convert(Vector{Float64}, collect(probs)), 1.0 - probsum)
        return new(elemsvector, probsvector)
    end
end

SamplerWeighted(elems, probs) = SamplerWeighted{eltype(elems)}(elems, probs)

Base.eltype(::Type{SamplerWeighted{T}}) where {T} = T
Base.eltype(::Type{SamplerUniform{T}}) where {T} = T

function Base.rand(rng::AbstractRNG, sp::SamplerWeighted)
    r = rand(rng)
    j = 1
    @inbounds cumulative_prob = sp.probs[j]
    while cumulative_prob < r
        j += 1
        @inbounds cumulative_prob += sp.probs[j]
    end
    return @inbounds sp.elems[j]
end

"""
    randseq{A::Alphabet, sp::Sampler, len::Integer)

Generate a LongSequence{A} of length `len` with elements drawn from
the given sampler.

# Example:
```
# Generate 1000-length RNA with 4% chance of N, 24% for A, C, G, or U
julia> sp = SamplerWeighted(rna"ACGUN", fill(0.24, 4))
julia> seq = randseq(RNAAlphabet{4}(), sp, 50)
50nt RNA Sequence:
CUNGGGCCCGGGNAAACGUGGUACACCCUGUUAAUAUCAACNNGCGCUNU
```
"""
function randseq(A::Alphabet, sp::Sampler, len::Integer)
    seq = LongSequence{typeof(A)}(len)
    @inbounds for i in 1:len
        letter = rand(sp)
        unsafe_setindex!(seq, letter, i)
    end
    return seq
end

"""
    randseq{A::Alphabet, len::Integer}

Generate a LongSequence{A} of length `len` from the specified alphabet, drawn
from the default distribution. User-defined alphabets should implement this
method to implement random LongSequence generation.

For RNA and DNA alphabets, the default distribution is uniform across A, C, G,
and T/U.
For AminoAcidAlphabet, it is uniform across the 20 proteogenic amino acids.
For a user-defined alphabet A, default is uniform across all elements of
`symbols(A)`.

# Example:
```
julia> seq = randseq(AminoAcidAlphabet(), 50)
50aa Amino Acid Sequence:
VFMHSIRMIRLMVHRSWKMHSARHVNFIRCQDKKWKSADGIYTDICKYSM
```
"""
function randseq(A::NucleicAcidAlphabet{4}, len::Integer)
    seq = LongSequence{typeof(A)}(len)
    randombits = zero(UInt64)
    # Generate random bases for each of the UInt64 in the .data field:
    for i in 1:seq_data_len(typeof(A), len)
        x = zero(UInt64)
        # 64 bits of randomness is enough for 2x16 4-bit nucleotides
        if isodd(i)
            randombits = rand(UInt64)
        end
        # Fill each of the 16 bases in the UInt64
        for nucleotide in 1:16
            # This is the position of the set bit of the four encoding the base
            bitposition = (randombits & 3) + 1
            x <<= bitposition
            x |= UInt64(1)
            x <<= 4 - bitposition
            randombits >>>= 2
        end
        @inbounds seq.data[i] = x
    end
    return seq
end

function randseq(::AminoAcidAlphabet, len::Integer)
    return randseq(AminoAcidAlphabet(), DefaultAASampler, len)
end

# Generic fallback method.
function randseq(A::Alphabet, len::Integer)
    if length(A) == 1 << bits_per_symbol(A)
        return randseq_allbits(A, len)
    else
        letters = symbols(A)
        sampler = SamplerUniform{eltype(A)}(letters)
        return randseq(A, sampler, len)
    end
end

function randseq_allbits(A::Alphabet, len::Integer)
    vector = rand(UInt64, seq_data_len(typeof(A), len))
    return LongSequence{typeof(A)}(vector, 1:len, false)
end

"""
    randdnaseq(len::Integer)

Generate a random LongSequence{DNAAlphabet{4}} sequence of length `len`,
with bases sampled uniformly from [A, C, G, T]
"""
randdnaseq(len::Integer) = randseq(DNAAlphabet{4}(), len)

"""
    randrnaseq(len::Integer)

Generate a random LongSequence{RNAAlphabet{4}} sequence of length `len`,
with bases sampled uniformly from [A, C, G, U]
"""
randrnaseq(len::Integer) = randseq(RNAAlphabet{4}(), len)

"""
    randaaseq(len::Integer)

Generate a random LongSequence{AminoAcidAlphabet} sequence of length `len`,
with amino acids sampled uniformly from the 20 standard amino acids.
"""
randaaseq(len::Integer) = randseq(AminoAcidAlphabet(), len)
