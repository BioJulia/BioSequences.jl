###
### Random Sequence Generator
###
###
### Random LongSequence generator.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

import Random: Sampler, rand!, GLOBAL_RNG

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
        elemsvector = convert(Vector{T}, vec(collect(elems)))
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
        elemsvector = convert(Vector{T}, vec(collect(elems)))
        probsvector = convert(Vector{Float64}, vec(collect(probs)))
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
    probs = sp.probs
    @inbounds cumulative_prob = probs[j]
    while cumulative_prob < r
        j += 1
        @inbounds cumulative_prob += probs[j]
    end
    return @inbounds sp.elems[j]
end
Base.rand(sp::SamplerWeighted) = rand(GLOBAL_RNG, sp)

###################### Generic longsequence methods ############################
# If no RNG is passed, use the global one
Random.rand!(seq::LongSequence) = rand!(GLOBAL_RNG, seq)
Random.rand!(seq::LongSequence, sp::Sampler) = rand!(GLOBAL_RNG, seq, sp)
randseq(A::Alphabet, len::Integer) = randseq(GLOBAL_RNG, A, len)
randseq(A::Alphabet, sp::Sampler, len::Integer) = randseq(GLOBAL_RNG, A, sp, len)
randdnaseq(len::Integer) = randdnaseq(GLOBAL_RNG, len)
randrnaseq(len::Integer) = randrnaseq(GLOBAL_RNG, len)
randaaseq(len::Integer) = randaaseq(GLOBAL_RNG, len)

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
function randseq(rng::AbstractRNG, A::Alphabet, len::Integer)
    rand!(rng, LongSequence{typeof(A)}(undef, len))
end

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
    rand!(rng, LongSequence{typeof(A)}(undef, len), sp)
end

# The generic method dispatches to `iscomplete`, since then we don't need
# to instantiate a sampler, and can use random bitpatterns
Random.rand!(rng::AbstractRNG, seq::LongSequence{A}) where A = rand!(rng, iscomplete(A()), seq)

####################### Implementations #######################################

# If given a sampler explicitly, just draw from that
function Random.rand!(rng::AbstractRNG, seq::LongSequence, sp::Sampler)
    @inbounds for i in eachindex(seq)
        letter = rand(rng, sp)
        seq[i] = letter
    end
    return seq
end

# 4-bit nucleotides's default distribution are equal among
# the non-ambiguous ones
function Random.rand!(rng::AbstractRNG, seq::LongSequence{<:NucleicAcidAlphabet{4}})
    data = seq.data
    rand!(rng, data)
    @inbounds for i in eachindex(data)
        nuc = 0x1111111111111111
        mask = data[i]
        nuc = ((nuc & mask) << 1) | (nuc & ~mask)
        mask >>>= 1
        nuc = ((nuc & mask) << 2) | (nuc & ~mask)
        data[i] = nuc
    end
    return seq
end

# Special AA implementation: Do not create the AA sampler, use the one
# that's already included.
Random.rand!(rng::AbstractRNG, seq::LongAA) = rand!(rng, seq, DefaultAASampler)

# All bitpatterns are valid - we use built-in RNG on the data vector.
function Random.rand!(rng::AbstractRNG, ::Val{true}, seq::LongSequence)
    rand!(rng, seq.data)
    seq
end

# Not all bitpatterns are valid - default to using a SamplerUniform
function Random.rand!(rng::AbstractRNG, ::Val{false}, seq::LongSequence)
    A = Alphabet(seq)
    letters = symbols(A)
    sampler = SamplerUniform{eltype(A)}(letters)
    rand!(rng, seq, sampler)
end

############################ Aliases for convenience ########################
"""
    randdnaseq([rng::AbstractRNG], len::Integer)

Generate a random LongSequence{DNAAlphabet{4}} sequence of length `len`,
with bases sampled uniformly from [A, C, G, T]
"""
randdnaseq(rng::AbstractRNG, len::Integer) = randseq(rng, DNAAlphabet{4}(), len)

"""
    randrnaseq([rng::AbstractRNG], len::Integer)

Generate a random LongSequence{RNAAlphabet{4}} sequence of length `len`,
with bases sampled uniformly from [A, C, G, U]
"""
randrnaseq(rng::AbstractRNG, len::Integer) = randseq(rng, RNAAlphabet{4}(), len)

"""
    randaaseq([rng::AbstractRNG], len::Integer)

Generate a random LongSequence{AminoAcidAlphabet} sequence of length `len`,
with amino acids sampled uniformly from the 20 standard amino acids.
"""
randaaseq(rng::AbstractRNG, len::Integer) = randseq(rng, AminoAcidAlphabet(), len)
