### BioSequences.jl
###
### A julia package for the representation and manipulation of biological sequences.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE

module BioSequences

export
    ###
    ### Symbols
    ###

    # Types & aliases
    NucleicAcid,
    DNA,
    RNA,
    DNA_A,
    DNA_C,
    DNA_G,
    DNA_T,
    DNA_M,
    DNA_R,
    DNA_W,
    DNA_S,
    DNA_Y,
    DNA_K,
    DNA_V,
    DNA_H,
    DNA_D,
    DNA_B,
    DNA_N,
    DNA_Gap,
    ACGT,
    ACGTN,
    RNA_A,
    RNA_C,
    RNA_G,
    RNA_U,
    RNA_M,
    RNA_R,
    RNA_W,
    RNA_S,
    RNA_Y,
    RNA_K,
    RNA_V,
    RNA_H,
    RNA_D,
    RNA_B,
    RNA_N,
    RNA_Gap,
    ACGU,
    ACGUN,
    AminoAcid,
    AA_A,
    AA_R,
    AA_N,
    AA_D,
    AA_C,
    AA_Q,
    AA_E,
    AA_G,
    AA_H,
    AA_I,
    AA_L,
    AA_K,
    AA_M,
    AA_F,
    AA_P,
    AA_S,
    AA_T,
    AA_W,
    AA_Y,
    AA_V,
    AA_O,
    AA_U,
    AA_B,
    AA_J,
    AA_Z,
    AA_X,
    AA_Term,
    AA_Gap,

    # Predicates
    isGC,
    iscompatible,
    isambiguous,
    iscertain,
    isgap,
    ispurine,
    ispyrimidine,

    ###
    ### Alphabets
    ###

    # Types & aliases
    Alphabet,
    NucleicAcidAlphabet,
    DNAAlphabet,
    RNAAlphabet,
    AminoAcidAlphabet,

    ###
    ### BioSequences
    ###

    join!,

    # Type & aliases
    BioSequence,
    NucSeq,
    AASeq,

    # Predicates
    ispalindromic,
    hasambiguity,
    isrepetitive,
    iscanonical,

    # Transformations
    canonical,
    canonical!,
    complement,
    complement!,
    reverse_complement,
    reverse_complement!,
    ungap,
    ungap!,
    join!,

    ###
    ### LongSequence
    ###

    # Type & aliases
    LongSequence,
    LongNuc,
    LongDNA,
    LongRNA,
    LongAA,
    LongSubSeq,

    # Random
    SamplerUniform,
    SamplerWeighted,
    randseq,
    randdnaseq,
    randrnaseq,
    randaaseq,

    ###
    ### Sequence literals
    ###

    @dna_str,
    @rna_str,
    @aa_str,

    @biore_str,
    @prosite_str,

    BioRegex,
    BioRegexMatch,
    matched,
    captured,
    alphabet, # TODO: Resolve the use of alphabet - it's from BioSymbols.jl
    symbols,
    gap,
    mismatches,
    matches,
    n_ambiguous,
    n_gaps,
    n_certain,
    gc_content,
    translate!,
    translate,
    ncbi_trans_table,

    # Search
    ExactSearchQuery,
    ApproximateSearchQuery,
    PFM,
    PWM,
    PWMSearchQuery,
    maxscore,
    scoreat,
    seqmatrix,
    majorityvote

using BioSymbols
import Twiddle: enumerate_nibbles,
    count_0000_nibbles,
    count_1111_nibbles,
    count_nonzero_nibbles,
    count_00_bitpairs,
    count_01_bitpairs,
    count_10_bitpairs,
    count_11_bitpairs,
    count_nonzero_bitpairs,
    repeatpattern
using Random

BioSymbols.gap(::Type{Char}) = '-'

include("alphabet.jl")

# Load the bit-twiddling internals that optimised BioSequences methods depend on.
include("bit-manipulation/bit-manipulation.jl")

# The generic, abstract BioSequence type
include("biosequence/biosequence.jl")

# The definition of the LongSequence concrete type, and its method overloads...

include("longsequences/longsequence.jl")
include("longsequences/hash.jl")
include("longsequences/randseq.jl")

include("geneticcode.jl")


# Pattern searching in sequences...
include("search/ExactSearchQuery.jl")
include("search/ApproxSearchQuery.jl")
include("search/re.jl")
include("search/pwm.jl")



struct Search{Q,I}
    query::Q
    itr::I
    overlap::Bool
end

const DEFAULT_OVERLAP = true

search(query, itr; overlap = DEFAULT_OVERLAP) = Search(query, itr, overlap)

function Base.iterate(itr::Search, state=firstindex(itr.itr))
    val = findnext(itr.query, itr.itr, state)
    val === nothing && return nothing
    state = itr.overlap ? first(val) + 1 : last(val) + 1
    return val, state
end

const HasRangeEltype = Union{<:ExactSearchQuery, <:ApproximateSearchQuery, <:Regex}

Base.eltype(::Type{<:Search{Q}}) where {Q<:HasRangeEltype} = UnitRange{Int}
Base.eltype(::Type{<:Search}) = Int
Base.IteratorSize(::Type{<:Search}) = Base.SizeUnknown()

function Base.findall(pat, seq::BioSequence; overlap::Bool = DEFAULT_OVERLAP)
    return collect(search(pat, seq; overlap))
end

function Base.findall(pat, seq::BioSequence, rng::UnitRange{Int}; overlap::Bool = DEFAULT_OVERLAP)
    v = view(seq, rng)
    itr = search(pat, v; overlap)
    return map(x->parentindices(v)[1][x], itr)
end

end  # module BioSequences
