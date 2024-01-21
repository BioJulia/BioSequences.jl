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
function LongSequence{A}(s::AbstractString) where {A <: Alphabet}
    return parse(LongSequence{A}, s)
end

function LongSequence{A}(
    src::Union{AbstractString,AbstractVector{UInt8}},
    part::AbstractUnitRange{<:Integer}=1:length(src)
) where {A<:Alphabet}
    len = length(part)
    seq = LongSequence{A}(undef, len)
    return copyto!(seq, 1, src, first(part), len)
end

Base.parse(::Type{T}, s::AbstractString) where {T <: LongSequence} = parse(T, String(s))

function Base.parse(T::Type{LongSequence{A}}, s::ASCIILike) where {A<:Alphabet}
    C = codetype(A())
    src = C isa AsciiAlphabet ? codeunits(s) : s
    n = _tryparse(T, s, C)
    if n isa Int
        throw_encode_error(A(), src, n)
    else
        n
    end
end

Base.tryparse(::Type{T}, s::AbstractString) where {T <: LongSequence} = tryparse(T, String(s))

function Base.tryparse(::Type{LongSequence{A}}, s::ASCIILike) where {A <: Alphabet}
    n = _tryparse(LongSequence{A}, s, codetype(A()))
    n isa Int ? nothing : n
end

function _tryparse(::Type{LongSequence{A}}, s::ASCIILike, ::AlphabetCode) where {A<:Alphabet}
    len = length(s)
    seq = LongSequence{A}(undef, len)
    # TODO!
    return copyto!(seq, 1, s, 1, len)
end

function _tryparse(::Type{LongSequence{A}}, s::ASCIILike, ::AsciiAlphabet) where {A<:Alphabet}
    seq = LongSequence{A}(undef, ncodeunits(s))
    try_encode_chunks!(seq, 1, codeunits(s), 1, ncodeunits(s))
end