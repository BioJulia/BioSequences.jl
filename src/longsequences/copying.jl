###
### Copying
###
###
### Copying methods for biological sequences.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

function Base.copyto!(dst::LongSequence{A}, src::LongSequence{A}) where {A}
    return copyto!(dst, 1, src, 1)
end

function Base.copyto!(dst::LongSequence{A}, doff::Integer,
                    src::LongSequence{A}, soff::Integer) where {A}
    return copyto!(dst, doff, src, soff, length(src) - soff + 1)
end

@inline function Base.copyto!(dst::LongSequence{A}, doff::Integer,
                    src::LongSequence{A}, soff::Integer, len::Integer) where {A}
    @boundscheck checkbounds(dst, doff:doff+len-1)
    @boundscheck checkbounds(src, soff:soff+len-1)

    if dst.shared || (dst === src && doff > soff)
        orphan!(dst, length(dst), true)
    end

    id = bitindex(dst, doff)
    is = bitindex(src, soff)
    #TODO: Resolve this use of bits_per_symbol
    rest = len * bits_per_symbol(A())
    dstdata = dst.data
    srcdata = src.data

    @inbounds while rest > 0
        # move `k` bits from `src` to `dst`
        x = dstdata[index(id)]
        y = srcdata[index(is)]
        if offset(id) < offset(is)
            y >>= offset(is) - offset(id)
            k = min(64 - offset(is), rest)
        else
            y <<= offset(id) - offset(is)
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

"""
    copy!(seq, src)

In-place copy content of `src` to biological sequence `seq`, resizing `seq` to fit.
`src` may be any sequence-like object that can be encoded to `seq`.

# Examples
```
julia> seq = copy!(dna"TAG", "AACGTM")
4nt DNA Sequence:
AACGTM

julia> copy!(seq, [0x61, 0x43, 0x54])
3nt DNA Sequence:
ACT
```
"""
function Base.copy!(dst::LongSequence{A}, src::LongSequence{A}) where {A <: Alphabet}
    len = length(src)
    resize!(dst, len, false)
    @inbounds copyto!(dst, 1, src, 1, len)
end

# Dispatch to alphabet type
function Base.copy!(seq::LongSequence{A}, src) where {A<:Alphabet}
   return copy!(seq, src, codetype(A()))
end

# Fast path for String + ASCII
function Base.copy!(seq::LongSequence, src::String, ::AsciiAlphabet)
   v = unsafe_wrap(Vector{UInt8}, src)
   resize!(seq, length(v), false)
   return encode_chunks!(seq, 1, v, 1, length(v))
end

# Generic method, cache len 'cause may be O(n) to calculate
function Base.copy!(seq::LongSequence, src, ::AlphabetCode)
   len = length(src)
   resize!(seq, len, false)
   return encode_copy!(seq, 1, src, 1, len)
end

# Actually, users don't need to create a copy of a sequence.
function Base.copy(seq::LongSequence)
    if seq.shared
        newseq = typeof(seq)(seq.data, seq.part, true)
        orphan!(newseq, length(seq), true)
    else
        newseq = typeof(seq)(copy(seq.data), seq.part, false)
    end
    return newseq
end

###
### Encoded copying
###

function encode_copy!(dst::LongSequence{A},
                      src::Union{AbstractVector,AbstractString}) where {A  <: Alphabet}
    return encode_copy!(dst, 1, src, 1)
end

function encode_copy!(dst::LongSequence{A},
                      doff::Integer,
                      src::Union{AbstractVector,AbstractString},
                      soff::Integer) where {A <: Alphabet}
    return encode_copy!(dst, doff, src, soff, codetype(A()))
end

function encode_copy!(dst::LongSequence{A},
                      doff::Integer,
                      src::Union{AbstractVector,AbstractString},
                      soff::Integer,
                      N::Integer) where {A <: Alphabet}
    return encode_copy!(dst, doff, src, soff, N, codetype(A()))
end


function encode_copy!(dst::LongSequence{A},
                      doff::Integer,
                      src::String,
                      soff::Integer, C::AsciiAlphabet) where {A <: Alphabet}
    v = unsafe_wrap(Vector{UInt8}, src)
    return encode_copy!(dst, doff, v, soff, length(v) - soff + 1, C)
end

function encode_copy!(dst::LongSequence{A},
                      doff::Integer,
                      src::Union{AbstractVector,AbstractString},
                      soff::Integer, C::AlphabetCode) where {A <: Alphabet}
    return encode_copy!(dst, doff, src, soff, length(src) - soff + 1, C)
end

function encode_copy!(dst::LongSequence{A},
                      doff::Integer,
                      src::Union{AbstractVector,AbstractString},
                      soff::Integer,
                      len::Integer, ::AlphabetCode) where {A}
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

@noinline function throw_encode_error(A::Alphabet, src::AbstractArray, soff::Integer)
    for i in 1:div(64, bits_per_symbol(A))
        sym = src[soff+i-1]
        stringbyte(A, sym) & 0x80 == 0x80 && error("Cannot encode $sym to $A")
    end
end

@inline function encode_chunk(A::Alphabet, src, soff::Int, N::Int)
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

# Use this for AsiiAlphabet alphabets only, internal use only, no boundschecks
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

@noinline function throw_enc_indexerr(N::Integer, len::Integer, soff::Integer)
    throw(ArgumentError("source of length $len does not contain $N elements from $soff"))
end

function encode_copy!(dst::LongSequence{A}, doff::Integer, src::AbstractVector{UInt8},
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
