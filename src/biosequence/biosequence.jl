###
### An abstract biological sequence type.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

abstract type BioSequence{A<:Alphabet} end

# Aliases and shorthands for describing subsets of the BioSequence type...
const NucleotideSeq = BioSequence{<:NucleicAcidAlphabet}
const AminoAcidSeq = BioSequence{AminoAcidAlphabet}

###
### Required traits and methods
###

# Base.length must be defined for every T<:BioSequence.

# As must the following...

"Get the length of a biological sequence."
@inline function Base.length(seq::BioSequence)
    error(
        string(
            "Base.length has not been defined for BioSequence type: ",
            typeof(seq),
            ". Any compatible concrete BioSequence subtype must have this method implemented."
        )
    )
end

"""
Return the data member of `seq` that stores the encoded sequence data.
"""
@inline function encoded_data(seq::BioSequence)
    error(
        string(
            "encoded_data has not been defined for BioSequence type: ",
            typeof(seq),
            ". Any compatible concrete BioSequence subtype must have this method implemented."
        )
    )
end

###
### Provided traits and methods
###

# These traits and methods are defined automatically for any subtype of BioSequence{A}.
# They may be overloaded for your concrete BioSequence sub-type if it is nessecery.

"Get the vector of bits storing a sequences packed encoded elements."
@inline encoded_data_type(seq::BioSequence) = typeof(encoded_data(seq))

"Get the element type of the vector of bits storing a sequences packed encoded elements."
@inline encoded_data_eltype(seq::BioSequence) = eltype(encoded_data_type(seq))

"""
Return the `Alpahbet` defining the possible biological symbols
and their encoding for a given biological sequence.
"""
@inline function Alphabet(::Type{<:BioSequence{A}}) where {A <: Alphabet}
    return A()
end

"Return the `Alpahbet` type that defines the biological symbols allowed for `seq`."
@inline function Alphabet(seq::BioSequence)
    return Alphabet(typeof(seq))
end

BioSymbols.alphabet(::Type{BioSequence{A}}) where {A<:Alphabet} = alphabet(A)

BitsPerSymbol(seq::BioSequence) = BitsPerSymbol(Alphabet(seq))

"Get the number of bits each symbol packed into a BioSequence uses, as an integer value."
bits_per_symbol(seq::BioSequence) = bits_per_symbol(Alphabet(seq))

# The generic functions for any BioSequence...
include("indexing.jl")
include("conversion.jl")
include("predicates.jl")
include("find.jl")
include("printing.jl")
include("transformations.jl")
include("counting.jl")