# site_counting.jl
# ================
#
# Site counting framework for biological sequences.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md


# Site types
# ----------

abstract type Site end
abstract type Position <: Site end

"""
A `Certain` site describes a site where both of two aligned sites are not an
ambiguity symbol or a gap.
"""
struct Certain <: Position end

"""
An `Gap` site describes a site where either of two aligned sites are a
gap symbol '-'.
"""
struct Gap <: Position end

"""
An `Ambiguous` site describes a site where either of two aligned sites are an
ambiguity symbol.
"""
struct Ambiguous <: Position end

"""
A `Match` site describes a site where two aligned nucleotides are the
same biological symbol.
"""
struct Match <: Position end

"""
A `Mismatch` site describes a site where two aligned nucleotides are not the
same biological symbol.
"""
struct Mismatch <: Position end


# Bitparallel counting
# --------------------
include("count_sites_bitpar.jl")

# Now we overload some internal methods for each S <: Site.
# So as different types of site can be counted in bit-parallel
# manner.

for A in (DNAAlphabet, RNAAlphabet)
    @eval begin

        # Gaps
        @inline function bp_chunk_count(::Type{Gap}, ::Type{$A{4}}, x::UInt64)
            return count_zero_nibbles(x)
        end

        @inline function bp_chunk_count(::Type{Gap}, ::Type{$A{4}}, a::UInt64, b::UInt64)
            # Count the gaps in a, count the gaps in b, subtract the number of shared gaps.
            return count_zero_nibbles(a) + count_zero_nibbles(b) - count_zero_nibbles(a | b)
        end

        # Certain
        @inline function bp_chunk_count(::Type{Certain}, ::Type{$A{4}}, x::UInt64)
            x = enumerate_nibbles(x)
            x = x ⊻ 0x1111111111111111
            return count_zero_nibbles(x)
        end

        @inline function bp_chunk_count(::Type{Certain}, ::Type{$A{4}}, a::UInt64, b::UInt64)
            x = enumerate_nibbles(a) ⊻ 0x1111111111111111
            y = enumerate_nibbles(b) ⊻ 0x1111111111111111
            return count_zero_nibbles(x | y)
        end

        # Ambiguous
        @inline function bp_chunk_count(::Type{Ambiguous}, ::Type{$A{4}}, x::UInt64)
            return count_nonzero_nibbles(enumerate_nibbles(x) & 0xEEEEEEEEEEEEEEEE)
        end

        @inline function bp_chunk_count(::Type{Ambiguous}, ::Type{$A{4}}, a::UInt64, b::UInt64)
            return count_nonzero_nibbles((enumerate_nibbles(a) | enumerate_nibbles(b)) & 0xEEEEEEEEEEEEEEEE)
        end

        # Match
        @inline function bp_chunk_count(::Type{Match}, ::Type{$A{4}}, a::UInt64, b::UInt64)
            return count_zero_nibbles(a ⊻ b)
        end

        @inline function bp_chunk_count(::Type{Match}, ::Type{$A{2}}, a::UInt64, b::UInt64)
            return count_zero_bitpairs(a ⊻ b)
        end

        # Mismatch
        @inline function bp_chunk_count(::Type{Mismatch}, ::Type{$A{4}}, a::UInt64, b::UInt64)
            return count_nonzero_nibbles(a ⊻ b)
        end

        @inline function bp_chunk_count(::Type{Mismatch}, ::Type{$A{2}}, a::UInt64, b::UInt64)
            return count_nonzero_bitpairs(a ⊻ b)
        end
    end
end

for s in (Match, Gap)
    @eval bp_correct_emptyspace(::Type{$s}, ::Type{A}) where {A<:NucleicAcidAlphabet} = true
end


# Specific Base.count methods
# ---------------------------

# In a general case, go to the bit-parallel counting method.
@inline Base.count(::Type{P}, a::MutableBioSequence{A}, b::MutableBioSequence{A}) where {P<:Position,A<:NucleicAcidAlphabet} = bitpar_counter(P, a, b)

# Some specific edge cases...
for A in (DNAAlphabet, RNAAlphabet)

    # Specific count methods for some edge cases regarding counting sites in
    # 2 bit encoded DNA and RNA sequences.
    seqtype = MutableBioSequence{A{2}}
    @eval begin
        @inline function Base.count(::Type{Certain}, a::$seqtype, b::$seqtype)
            return min(length(a), length(b))
        end
        @inline Base.count(::Type{Gap}, a::$seqtype, b::$seqtype) = 0
        @inline Base.count(::Type{Ambiguous}, a::$seqtype, b::$seqtype) = 0
    end
end

# Specific Base.count sliding window methods
# ------------------------------------------

function Base.count(::Type{P}, a::MutableBioSequence{A}, b::MutableBioSequence{A}, width::Int, step::Int) where {P<:Position,A<:NucleicAcidAlphabet}
    len = min(length(a), length(b))
    ritr = StepRange(width, step, len)
    width -= 1
    results = Vector{IntervalValue{Int,bp_counter_type(P, A)}}(undef, length(ritr))
    r = 1
    @inbounds for i in ritr
        idx = (i - width):i
        results[r] = IntervalValue(first(idx), last(idx), count(P, a[idx], b[idx]))
        r += 1
    end
    return results
end

# Specific count_pairwise methods
# -------------------------------

@inline diag_val(::Type{Int}) = Int(0)

function count_pairwise(::Type{P}, seqs::Vararg{MutableBioSequence{A},N}) where {P<:Position,A<:NucleicAcidAlphabet,N}
    @assert N >= 2 "At least two sequences are required."
    counts = Matrix{bp_counter_type(P, A)}(undef, (N, N))
    for i in 1:N
        counts[i,i] = diag_val(eltype(counts))
        for j in (i+1):N
            counts[i,j] = counts[j,i] = count(P, seqs[i], seqs[j])
        end
    end
    return counts
end


# General, non-specific Base.count methods
# ----------------------------------------

"""
    Base.count{S<:Site}(::Type{S}, a::BioSequence, b::BioSequence)

Count the number of sites of type `S`, between two sequences a and b.
"""
@inline function Base.count(::Type{S}, a::MutableBioSequence, b::MutableBioSequence) where {S<:Site}
    seqs = promote(a, b)
    return count(S, seqs...)
end

"""
    Base.count{S<:Site}(::Type{S}, a::BioSequence, b::BioSequence, width::Int, step::Int)

Count the number of sites of type `S`, between two sequences a and b, using a
sliding window of a given `width`, and `step`.
The `width` and `step` must be provided as number of base pairs.
"""
function Base.count(::Type{S}, a::MutableBioSequence, b::MutableBioSequence, width::Int, step::Int) where {S<:Site}
    seqs = promote(a, b)
    return count(S, seqs..., width, step)
end

"""
    count_pairwise{S<:Site,N}(::Type{S}, seqs::Vararg{BioSequence,N})

Count the number of sites of type `S`, between each possible pair of the sequences
provided (`seqs`).

This returns a symmetric matrix contining the output for each pairwise
operation.
"""
function count_pairwise(::Type{S}, seqs::Vararg{MutableBioSequence,N}) where {S<:Site,N}
    seqs = promote(seqs...)
    return count_pairwise(S, seqs...)
end
