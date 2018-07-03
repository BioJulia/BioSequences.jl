```@meta
CurrentModule = BioSequences
DocTestSetup = quote
    using BioSequences
end
```
# [The abstract `BioSequence` type: Developer Information](@id biosequence_dev)

docs@```
BioSequence
```

### On concrete biological sequence types...

Sequences in BioSequences.jl are more strictly typed than in many other
libraries;
Elements in a sequence are typed as biological symbol instead of character or
byte. They are special-purpose types rather than simply strings and hence
offer additional functionality that naive string types don't have.

Though this strictness sacrifices some convenience, it also means you can
always rely on a DNA sequence type to store DNA and nothing but DNA, without
having to check, or deal with lowercase versus uppercase and so on.

Strict separation of sequence types also means we are free to choose the most
efficient representation. DNA and RNA sequences are encoded using either four
bits per base (which is the default), or two bits per base. This makes them
memory efficient and allows us to speed up many common operations and
transformations, like nucleotide composition, reverse complement, and *k*-mer
enumeration.


# [Alphabet traits](@id alphabets)

docs@```
Alphabet
```


## Defining a new alphabet

The alphabet type parameter `A` of `BioSequence{A}` enables a user to extend
functionality of `BioSequence` with minimum effort. As an example, definition of
a new alphabet type representing a sequence of boolean values is shown below:

```jldoctest
julia> immutable BoolAlphabet <: Alphabet end

julia> BioSequences.bitsof(::Type{BoolAlphabet}) = 1

julia> BioSequences.eltype(::Type{BoolAlphabet}) = Bool

julia> BioSequences.alphabet(::Type{BoolAlphabet}) = false:true

julia> function BioSequences.encode(::Type{BoolAlphabet}, x::Bool)
           return UInt64(ifelse(x, 0x01, 0x00))
       end

julia> function BioSequences.decode(::Type{BoolAlphabet}, x::UInt64)
           if x > 0x01
               throw(BioSequences.DecodeError(BoolAlphabet, x))
           end
           return ifelse(x == 0x00, false, true)
       end

```
