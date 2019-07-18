###
### Random Sequence Generator
###
###
### Random sequence generator.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

import Random: Sampler, GLOBAL_RNG

"""
    SamplerUniform{T}

Uniform sampler of type T. Instantiate with a collection of eltype T containing
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
        if length(elemsvector) < 1
            throw(ArgumentError("elements collection must be non-empty"))
        end
        return new(elemsvector)
    end
end

SamplerUniform(elems) = SamplerUniform{eltype(elems)}(elems)
Base.eltype(::Type{SamplerUniform{T}}) where {T} = T
Base.rand(rng::AbstractRNG, sp::SamplerUniform) = rand(rng, sp.elems)
Base.rand(sp::SamplerUniform) = rand(GLOBAL_RNG, sp.elems)
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
        elemsvector = convert(Vector{T}, collect(elems))
        probsvector = convert(Vector{Float64}, collect(probs))
        if !isempty(probsvector)
            probsum = sum(probsvector)
            if probsum > 1.0
                throw(ArgumentError("sum of probabilties cannot exceed 1.0"))
            elseif minimum(probsvector) < 0.0
                throw(ArgumentError("probabilities must be non-negative"))
            end
        else
            probsum = 0.0
        end
        if length(elemsvector) != length(probsvector) + 1
            throw(ArgumentError("length of elems must be length of probs + 1"))
        end
        # Even with float weirdness,  we can guarantee x + (1.0 - x) == 1.0,
        # when 0 ≤ x ≤ 1, as there's no exponent, and it works like int addition
        return new(elemsvector, push!(probsvector, 1.0 - probsum))
    end
end

SamplerWeighted(elems, probs) = SamplerWeighted{eltype(elems)}(elems, probs)
Base.eltype(::Type{SamplerWeighted{T}}) where {T} = T

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
Base.rand(sp::SamplerWeighted) = rand(GLOBAL_RNG, sp)

"""
    randseq([rng::AbstractRNG], A::Alphabet, sp::Sampler, len::Integer)

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
function randseq(rng::AbstractRNG, A::Alphabet, sp::Sampler, len::Integer)
    seq = LongSequence{typeof(A)}(len)
    @inbounds for i in 1:len
        letter = rand(rng, sp)
        unsafe_setindex!(seq, letter, i)
    end
    return seq
end
randseq(A::Alphabet, sp::Sampler, len::Integer) = randseq(GLOBAL_RNG, A, sp, len)

"""
    randseq([rng::AbstractRNG], A::Alphabet, len::Integer)

Generate a LongSequence{A} of length `len` from the specified alphabet, drawn
from the default distribution. User-defined alphabets should implement this
method to implement random LongSequence generation.

For RNA and DNA alphabets, the default distribution is uniform across A, C, G,
and T/U.
For AminoAcidAlphabet, it is uniform across the 20 standard amino acids.
For a user-defined alphabet A, default is uniform across all elements of
`symbols(A)`.

# Example:
```
julia> seq = randseq(AminoAcidAlphabet(), 50)
50aa Amino Acid Sequence:
VFMHSIRMIRLMVHRSWKMHSARHVNFIRCQDKKWKSADGIYTDICKYSM
```
"""
function randseq(rng::AbstractRNG, A::NucleicAcidAlphabet{4}, len::Integer)
    seq = LongSequence{typeof(A)}(len)
    randombits = zero(UInt64)
    # Generate random bases for each of the UInt64 in the .data field:
    for i in 1:seq_data_len(typeof(A), len)
        x = zero(UInt64)
        # 64 bits of randomness is enough for 2x16 4-bit nucleotides
        if isodd(i)
            randombits = rand(rng, UInt64)
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
randseq(A::NucleicAcidAlphabet{4}, len::Integer) = randseq(GLOBAL_RNG, A, len)

function randseq(rng::AbstractRNG, ::AminoAcidAlphabet, len::Integer)
    return randseq(rng, AminoAcidAlphabet(), DefaultAASampler, len)
end
randseq(A::AminoAcidAlphabet, len::Integer) = randseq(GLOBAL_RNG, A, len)

# Generic fallback method.
function randseq(rng::AbstractRNG, A::Alphabet, len::Integer)
    if length(A) == 1 << bits_per_symbol(A)
        return randseq_allbits(rng, A, len)
    else
        letters = symbols(A)
        sampler = SamplerUniform{eltype(A)}(letters)
        return randseq(rng, A, sampler, len)
    end
end
randseq(A::Alphabet, len::Integer) = randseq(GLOBAL_RNG, A, len)

function randseq_allbits(rng::AbstractRNG, A::Alphabet, len::Integer)
    vector = rand(rng, UInt64, seq_data_len(typeof(A), len))
    return LongSequence{typeof(A)}(vector, 1:len, false)
end

"""
    randdnaseq([rng::AbstractRNG], len::Integer)

Generate a random LongSequence{DNAAlphabet{4}} sequence of length `len`,
with bases sampled uniformly from [A, C, G, T]
"""
randdnaseq(rng::AbstractRNG, len::Integer) = randseq(rng, DNAAlphabet{4}(), len)
randdnaseq(len::Integer) = randseq(GLOBAL_RNG, DNAAlphabet{4}(), len)

"""
    randrnaseq([rng::AbstractRNG], len::Integer)

Generate a random LongSequence{RNAAlphabet{4}} sequence of length `len`,
with bases sampled uniformly from [A, C, G, U]
"""
randrnaseq(rng::AbstractRNG, len::Integer) = randseq(rng, RNAAlphabet{4}(), len)
randrnaseq(len::Integer) = randseq(GLOBAL_RNG, RNAAlphabet{4}(), len)

"""
    randaaseq([rng::AbstractRNG], len::Integer)

Generate a random LongSequence{AminoAcidAlphabet} sequence of length `len`,
with amino acids sampled uniformly from the 20 standard amino acids.
"""
randaaseq(rng::AbstractRNG, len::Integer) = randseq(rng, AminoAcidAlphabet(), len)
randaaseq(len::Integer) = randseq(GLOBAL_RNG, AminoAcidAlphabet(), len)
