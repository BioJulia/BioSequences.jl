###
### Copying
###
###
### Copying methods for biological sequences.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

# TODO: Add generic emethods for other bioseqs like mers and refseqs
##########
"""
    copy!(dst::LongSequence, src::BioSequence)

In-place copy content of `src` to `dst`, resizing `dst` to fit.
The alphabets of `src` and `dst` must be compatible.

# Examples
```
julia> seq = copy!(dna"TAG", dna"AACGTM")
6nt DNA Sequence:
AACGTM

julia> copy!(seq, rna"UUAG")
4nt DNA Sequence:
TTAG
```
"""
function Base.copy!(dst::LongSequence{A}, src::LongSequence{A}) where {A <: Alphabet}
    return _copy!(dst, src)
end

function Base.copy!(dst::LongSequence{<:NucleicAcidAlphabet{N}},
                    src::LongSequence{<:NucleicAcidAlphabet{N}}) where N
    return _copy!(dst, src)
end

function _copy!(dst::LongSequence, src::LongSequence)
    orphan!(dst)
    # Most efficient
    if !src.shared
        resize!(dst.data, length(src.data))
        copyto!(dst.data, src.data)
        dst.part = src.part
    # Less efficient
    else
        resize!(dst, length(src))
        copyto!(dst, src)
    end
    return dst
end

"""
    copyto!(dst::LongSequence, src::BioSequence)

Equivalent to `copyto!(dst, 1, src, 1, length(src))`
"""
function Base.copyto!(dst::LongSequence{A}, src::LongSequence{A}) where {A <: Alphabet}
    return copyto!(dst, 1, src, 1, length(src))
end

function Base.copyto!(dst::LongSequence{<:NucleicAcidAlphabet{N}},
                      src::LongSequence{<:NucleicAcidAlphabet{N}}) where N
    return copyto!(dst, 1, src, 1, length(src))
end

"""
    copyto!(dst::LongSequence, soff, src::BioSequence, doff, N)

In-place copy `N` elements from `src` starting at `soff` to `dst`, starting at `doff`.
The length of `dst` must be greater than or equal to `N + doff - 1`.
The first N elements of `dst` are overwritten,
the other elements are left untouched. The alphabets of `src` and `dst` must be compatible.

# Examples
```
julia> seq = copyto!(dna"AACGTM", 1, dna"TAG", 1, 3)
6nt DNA Sequence:
TAGGTM

julia> copyto!(seq, 2, rna"UUUU", 1, 4)
6nt DNA Sequence:
TTTTTM
```
"""
function Base.copyto!(dst::LongSequence{A}, doff::Integer,
                      src::LongSequence{A}, soff::Integer,
                      N::Integer) where {A <: Alphabet}
    return _copyto!(dst, doff, src, soff, N)
end

function Base.copyto!(dst::LongSequence{<:NucleicAcidAlphabet{B}}, doff::Integer,
                      src::LongSequence{<:NucleicAcidAlphabet{B}}, soff::Integer,
                      N::Integer) where B
    return _copyto!(dst, doff, src, soff, N)
end

function _copyto!(dst::LongSequence{A}, doff::Integer,
                  src::LongSequence, soff::Integer,
                  N::Integer) where {A <: Alphabet}
    @boundscheck checkbounds(dst, doff:doff+N-1)
    @boundscheck checkbounds(src, soff:soff+N-1)

    if dst.shared
        orphan!(dst)
    end
    # This prevents a sequence from destructively overwriting its own data
    if (dst === src) & (doff > soff)
        return _copyto!(dst, doff, copy(src[soff:soff+N-1]), 1, N)
    end

    id = bitindex(dst, doff)
    is = bitindex(src, soff)
    rest = N * bits_per_symbol(A())
    dstdata = dst.data
    srcdata = src.data

    @inbounds while rest > 0
        # move `k` bits from `src` to `dst`
        x = dstdata[index(id)]
        y = srcdata[index(is)]
        if offset(id) < offset(is)
            y >>= (offset(is) - offset(id)) & 63
            k = min(64 - offset(is), rest)
        else
            y <<= (offset(id) - offset(is)) & 63
            k = min(64 - offset(id), rest)
        end
        m = bitmask(k) << offset(id)
        dstdata[index(id)] = y & m | x & ~m

        id += k
        is += k
        rest -= k
    end

    return dst
end

function Base.copy(seq::LongSequence)
    if seq.shared
        newseq = typeof(seq)(seq.data, seq.part, true)
        orphan!(newseq, length(seq), true)
    else
        newseq = typeof(seq)(copy(seq.data), seq.part, false)
    end
    return newseq
end

#########
const SeqLike = Union{AbstractVector, AbstractString}
const ASCIILike = Union{String, SubString{String}}

"""
    copy!(dst::LongSequence, src)

In-place copy content of sequence-like object `src` to `dst`, resizing `dst` to fit.
The content of `src` must be able to be encoded to the alphabet of `dst`.

# Examples
```
julia> seq = copy!(dna"TAG", "AACGTM")
6nt DNA Sequence:
AACGTM

julia> copy!(seq, [0x61, 0x43, 0x54])
3nt DNA Sequence:
ACT
```
"""
function Base.copy!(dst::LongSequence{A}, src::SeqLike) where {A <: Alphabet}
    return copy!(dst, src, codetype(A()))
end

function Base.copy!(dst::LongSequence{<:Alphabet}, src::ASCIILike, C::AsciiAlphabet)
    v = GC.@preserve src unsafe_wrap(Vector{UInt8}, pointer(src), ncodeunits(src))
    return copy!(dst, v, C)
end

function Base.copy!(dst::LongSequence{<:Alphabet}, src::AbstractVector{UInt8}, ::AsciiAlphabet)
    resize!(dst, length(src))
    return encode_chunks!(dst, 1, src, 1, length(src))
end

function Base.copy!(dst::LongSequence{<:Alphabet}, src::SeqLike, ::AlphabetCode)
    len = length(src) # calculate only once
    resize!(dst, len)
    return copyto!(dst, 1, src, 1, len)
end


########

for (anum, atype) in enumerate((DNAAlphabet{4}, DNAAlphabet{2}, RNAAlphabet{4},
    RNAAlphabet{2}, AminoAcidAlphabet))
    tablename = Symbol("BYTE_TO_ALPHABET_CHAR" * string(anum))
    @eval begin
        alph = $(atype)()
        syms = symbols(alph)
        const $(tablename) = let
            bytes = fill(0x80, 256)
            for symbol in syms
                bytes[UInt8(Char(symbol)) + 1] = encode(alph, symbol)
                bytes[UInt8(lowercase(Char(symbol))) + 1] = encode(alph, symbol)
            end
            Tuple(bytes)
        end
        stringbyte(::$(atype), x::UInt8) = @inbounds $(tablename)[x + 1]
    end
end

# This is used to effectively scan an array of UInt8 for invalid bytes, when one is detected
@noinline function throw_encode_error(A::Alphabet, src::AbstractArray{UInt8}, soff::Integer)
    for i in 1:div(64, bits_per_symbol(A))
        sym = src[soff+i-1]
        stringbyte(A, sym) & 0x80 == 0x80 && error("Cannot encode $sym to $A")
    end
end

@inline function encode_chunk(A::Alphabet, src::AbstractArray{UInt8}, soff::Integer, N::Integer)
    chunk = zero(UInt64)
    check = 0x00
    @inbounds for i in 1:N
        enc = stringbyte(A, src[soff+i-1])
        check |= enc
        chunk |= UInt64(enc) << (bits_per_symbol(A) * (i-1))
    end
    check & 0x80 == 0x00 || throw_encode_error(A, src, soff)
    return chunk
end

# Use this for AsiiAlphabet alphabets only, internal use only, no boundschecks.
# This is preferential to `copyto!` if none of the sequence's original content
# needs to be kept, since this is faster.
function encode_chunks!(dst::LongSequence{A}, startindex::Integer, src::AbstractVector{UInt8},
                        soff::Integer, N::Integer) where {A <: Alphabet}
    chunks, rest = divrem(N, symbols_per_data_element(dst))
    @inbounds for i in startindex:startindex+chunks-1
        dst.data[i] = encode_chunk(A(), src, soff, symbols_per_data_element(dst))
        soff += symbols_per_data_element(dst)
    end
    @inbounds if !iszero(rest)
        dst.data[startindex+chunks] = encode_chunk(A(), src, soff, rest)
    end
    return dst
end

#########

# Two-argument method
"""
    copyto!(dst::LongSequence, src)

Equivalent to `copyto!(dst, 1, src, 1, length(src))`.
"""
function Base.copyto!(dst::LongSequence{A}, src::SeqLike) where {A <: Alphabet}
    return copyto!(dst, src, codetype(A()))
end

# Specialized method to avoid O(N) length call for string-like src
function Base.copyto!(dst::LongSequence{<:Alphabet}, src::ASCIILike, C::AsciiAlphabet)
    return copyto!(dst, 1, src, 1, ncodeunits(src), C)
end

function Base.copyto!(dst::LongSequence{<:Alphabet}, src::SeqLike, C::AlphabetCode)
    return copyto!(dst, 1, src, 1, length(src), C)
end

"""
    copyto!(dst::LongSequence, soff, src, doff, N)

In-place encode `N` elements from `src` starting at `soff` to `dst`, starting at `doff`.
The length of `dst` must be greater than or equal to `N + doff - 1`.
The first N elements of `dst` are overwritten,
the other elements are left untouched. The content of `src` must be able to be encoded to
the alphabet of `dst`.

# Examples
```
julia> seq = copyto!(dna"AACGTM", 1, "TAG", 1, 3)
6nt DNA Sequence:
TAGGTM

julia> copyto!(seq, 2, rna"UUUU", 1, 4)
6nt DNA Sequence:
TTTTTM
```
"""
# Dispatch to codetype
function Base.copyto!(dst::LongSequence{A}, doff::Integer,
                        src::SeqLike, soff::Integer, N::Integer) where {A <: Alphabet}
    return copyto!(dst, doff, src, soff, N, codetype(A()))
end

# For ASCII seq and src, convert to byte vector and dispatch using that
function Base.copyto!(dst::LongSequence{A}, doff::Integer,
                      src::ASCIILike, soff::Integer,
                      N::Integer, C::AsciiAlphabet) where {A <: Alphabet}
    v = GC.@preserve src unsafe_wrap(Vector{UInt8}, pointer(src), ncodeunits(src))
    return Base.copyto!(dst, doff, v, soff, N, C)
end

@noinline function throw_enc_indexerr(N::Integer, len::Integer, soff::Integer)
    throw(ArgumentError("source of length $len does not contain $N elements from $soff"))
end

# Generic method for copyto!, i.e. NOT ASCII input
function Base.copyto!(dst::LongSequence{A}, doff::Integer,
                      src::SeqLike, soff::Integer,
                      len::Integer, ::AlphabetCode) where {A <: Alphabet}
    if soff != 1 && isa(src, AbstractString) && !isascii(src)
        throw(ArgumentError("source offset ≠ 1 is not supported for non-ASCII string"))
    end

    checkbounds(dst, doff:doff+len-1)
    length(src) < soff + len - 1 && throw_enc_indexerr(len, length(src), soff)

    orphan!(dst)

    next = bitindex(dst, doff)
    stop = bitindex(dst, doff + len)

    i = soff
    while next < stop
        x = UInt64(0)
        j = index(next)
        while index(next) == j && next < stop
            char = src[i]
            i = nextind(src, i)
            x |= enc64(dst, convert(Char, char)) << offset(next)
            #TODO: Resolve this use of bits_per_symbol...
            next += bits_per_symbol(A())
        end
        dst.data[j] = x
    end
    return dst
end

# Special method possible for ASCII alphabet and UInt8 array
function Base.copyto!(dst::LongSequence{A}, doff::Integer, src::AbstractVector{UInt8},
                      soff::Integer, N::Integer, ::AsciiAlphabet) where {A<:Alphabet}
    checkbounds(dst, doff:doff+N-1)
    length(src) < soff + N - 1 && throw_enc_indexerr(N, length(src), soff)
    orphan!(dst)
    bitind = bitindex(dst, doff)
    remaining = N

    # Fill in first chunk. Since it may be only partially filled from both ends,
    # bit-tricks are harder and we do it the old-fashined way.
    # Maybe this can be optimized but eeehhh...
    @inbounds while (!iszero(offset(bitind)) & !iszero(remaining))
        dst[doff] = eltype(dst)(reinterpret(Char, (src[soff] % UInt32) << 24))
        doff += 1
        soff += 1
        remaining -= 1
        bitind += bits_per_symbol(A())
    end

    # Fill in middle
    n = remaining - rem(remaining, symbols_per_data_element(dst))
    encode_chunks!(dst, index(bitind), src, soff, n)
    remaining -= n
    soff += n
    bitind += n * bits_per_symbol(A())

    # Fill in last chunk
    @inbounds if !iszero(remaining)
        chunk = encode_chunk(A(), src, soff, remaining)
        before = dst.data[index(bitind)] & (typemax(UInt) << (remaining * bits_per_symbol(A())))
        dst.data[index(bitind)] = chunk | before
    end

    return dst
end
