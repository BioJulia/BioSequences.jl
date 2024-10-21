###
### Constructors
###
###
### Constructor methods for LongSequences.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

@inline seq_data_len(s::LongSequence{A}) where A = seq_data_len(A, length(s))

@inline function seq_data_len(::Type{A}, len::Integer)::Int where A <: Alphabet
    iszero(bits_per_symbol(A())) && return 0
    return cld(len % UInt, div(64, bits_per_symbol(A())) % UInt) % Int
end

function LongSequence{A}(::UndefInitializer, len::Integer) where {A<:Alphabet}
    if len < 0
        throw(ArgumentError("len must be non-negative"))
    end
    return LongSequence{A}(Memory{UInt64}(undef, seq_data_len(A, len)), UInt(len))
end

# Generic constructor
function LongSequence{A}(it) where {A <: Alphabet}
    len = length(it)
    data = Memory{UInt64}(undef, seq_data_len(A, len))
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

Base.empty(::Type{T}) where {T <: LongSequence} = T(Memory{UInt64}(), UInt(0))
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