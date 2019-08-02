### BioSequences.jl
###
### A julia package for the representation and manipulation of biological sequences.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

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
    CharAlphabet,
    
    ###
    ### BioSequences
    ###
    
    # Type & aliases
    BioSequence,
    NucleotideSeq,
    AminoAcidSeq,
    
    # Indexing
    unsafe_setindex!,
    
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
    
    # Iteration
    each,
    fwmer,
    bwmer,
    
    ###
    ### LongSequence
    ###
    
    # Type & aliases
    LongSequence,
    LongDNASeq,
    LongRNASeq,
    LongAminoAcidSeq,
    LongCharSeq,
    
    # Random
    SamplerUniform,
    SamplerWeighted,
    randseq,
    randdnaseq,
    randrnaseq,
    randaaseq,
    
    ###
    ### Mers
    ###
    
    # Type & aliases
    AbstractMer,
    Mer,
    DNAMer,
    RNAMer,
    DNAKmer,
    RNAKmer,
    
    BigMer,
    BigDNAMer,
    BigRNAMer,
    BigDNAKmer,
    BigRNAKmer,
    
    DNACodon,
    RNACodon,
    
    # Iteration
    neighbors,
    
    ###
    ### Sequence literals
    ###
    
    @dna_str,
    @rna_str,
    @aa_str,
    @char_str,
    @biore_str,
    @prosite_str,
    @mer_str,
    @bigmer_str,
    
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
    
    eachcanonical,
    
    ###
    ### Composition
    ###
    Composition,
    composition,
    NucleicAcidCounts,
    
    translate,
    ncbi_trans_table,
    
    
    # Search
    ExactSearchQuery,
    ApproximateSearchQuery,
    approxsearch,
    approxsearchindex,
    approxrsearch,
    approxrsearchindex,
    PFM,
    PWM,
    maxscore,
    scoreat,
    
    ReferenceSequence,
    
    ###
    ### Demultiplexing
    ###
    Demultiplexer,
    demultiplex,
    
    seqmatrix,
    majorityvote,
    MinHashSketch,
    minhash,
    Site,
    Certain,
    Ambiguous,
    Gap,
    Match,
    Mismatch,
    count_pairwise

using BioGenerics
using BioSymbols
import Combinatorics
import IndexableBitVectors
import Twiddle: enumerate_nibbles,
    nibble_mask,
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

# The definition of the Skipmer concrete type, and its method overloads...
include("mers/mer.jl")

# The definition of the ReferenceSequence concrete type, and its method overloads...
include("nmask.jl")
include("refseq/refseq.jl")

# The generic iterators for any BioSequence...
include("iterators/condition.jl")
include("iterators/eachmer.jl")
include("iterators/skipmerfactory.jl")

include("composition.jl")

include("geneticcode.jl")
include("demultiplexer.jl")

# Pattern searching in sequences...
include("search/exact.jl")
include("search/approx.jl")
include("search/re.jl")
include("search/pwm.jl")

include("minhash.jl")

end  # module BioSequences
