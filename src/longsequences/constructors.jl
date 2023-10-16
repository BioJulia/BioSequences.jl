###
### Constructors
###
###
### Constructor methods for LongSequences.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

@inline seq_data_len(s::LongSequence{A}) where A = seq_data_len(A, length(s))

@inline function seq_data_len(::Type{A}, len::Integer) where A <: Alphabet
    iszero(bits_per_symbol(A())) && return 0
    return cld(len, div(64, bits_per_symbol(A())))
end

function LongSequence{A}(::UndefInitializer, len::Integer) where {A<:Alphabet}
    if len < 0
        throw(ArgumentError("len must be non-negative"))
    end
    return LongSequence{A}(Vector{UInt64}(undef, seq_data_len(A, len)), UInt(len))
end

# Generic constructor
function LongSequence{A}(it) where {A <: Alphabet}
    len = length(it)
    data = Vector{UInt64}(undef, seq_data_len(A, len))
    bits = zero(UInt)
    bitind = bitindex(BitsPerSymbol(A()), encoded_data_eltype(LongSequence{A}), 1)
    @inbounds for x in it
        xT = convert(eltype(A), x)
        enc = encode(A(), xT)
        bits |= enc << offset(bitind)
        if iszero(offset(nextposition(bitind)))
            data[index(bitind)] = bits
            bits = zero(UInt64)
        end
        bitind = nextposition(bitind)
    end
    iszero(offset(bitind)) || (data[index(bitind)] = bits)
    LongSequence{A}(data, len % UInt)
end

Base.empty(::Type{T}) where {T <: LongSequence} = T(UInt[], UInt(0))
(::Type{T})() where {T <: LongSequence} = empty(T)

# Constructors from other sequences
# TODO: Remove this method, since the user can just slice
LongSequence(s::LongSequence, xs::AbstractUnitRange{<:Integer}) = s[xs]

function LongSequence(seq::BioSequence{A}) where {A <: Alphabet}
    return LongSequence{A}(seq)
end

LongSequence{A}(seq::LongSequence{A}) where {A <: Alphabet} = copy(seq)

function (::Type{T})(seq::LongSequence{<:NucleicAcidAlphabet{N}}) where
         {N, T<:LongSequence{<:NucleicAcidAlphabet{N}}}
    return T(copy(seq.data), seq.len)
end

# Constructors from strings
function LongSequence{A}(s::Union{String, SubString{String}}) where {A<:Alphabet}
    return LongSequence{A}(s, codetype(A()))
end

# Generic method for String/Substring.
function LongSequence{A}(s::Union{String, SubString{String}}, ::AlphabetCode) where {A<:Alphabet}
    len = length(s)
    seq = LongSequence{A}(undef, len)
    return copyto!(seq, 1, s, 1, len)
end

function LongSequence{A}(s::Union{String, SubString{String}}, ::AsciiAlphabet) where {A<:Alphabet}
    seq = LongSequence{A}(undef, ncodeunits(s))
    return encode_chunks!(seq, 1, codeunits(s), 1, ncodeunits(s))
end

function LongSequence{A}(
    src::Union{AbstractString,AbstractVector{UInt8}},
    part::AbstractUnitRange{<:Integer}=1:length(src)
) where {A<:Alphabet}
    len = length(part)
    seq = LongSequence{A}(undef, len)
    return copyto!(seq, 1, src, first(part), len)
end

Base.parse(::Type{LongSequence{A}}, seq::AbstractString) where A = LongSequence{A}(seq)

"""
    guess_parse(s::Union{AbstractString, AbstractVector{UInt8}}) -> LongSequence

Parse `s` into a `LongSequence` with an appropriate `Alphabet`, or throw an exception
if no alphabet matches.
See [`guess_alphabet`](@ref) for the available alphabets and the alphabet priority.

!!! warning
    The functions `guess_parse` and `guess_alphabet` are intended for use in interactive
    sessions, and are not suitable for use in packages or non-ephemeral work.
    They are type unstable, and their heuristics **are subject to change** in minor versions.

# Examples
```jldoctest
julia> guess_parse("QMKLPEEFW")
9aa Amino Acid Sequence:
QMKLPEEFW

julia> guess_parse("UAUGCUGUAGG")
11nt RNA Sequence:
UAUGCUGUAGG

julia> guess_parse("PKMW#3>>0;kL")
ERROR: cannot encode 0x23 in AminoAcidAlphabet
[...]
```
"""
function guess_parse(s::AbstractVector{UInt8})
    A = guess_alphabet(s)
    A isa Integer && throw(EncodeError(AminoAcidAlphabet(), s[A]))
    LongSequence{typeof(A)}(s)
end
guess_parse(s::Union{String, SubString{String}}) = guess_parse(codeunits(s))
guess_parse(s::AbstractString) = guess_parse(String(s))