# Copying
# =======
#
# Copying methods for biological sequences.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

function Base.copyto!(seq::MutableBioSequence{A}, doff::Integer,
                    src::Vector{UInt8},  soff::Integer, len::Integer) where {A}
    datalen = seq_data_len(A, len)
    return seq
end

function Base.copyto!(dst::MutableBioSequence{A}, src::MutableBioSequence{A}) where {A}
    return copyto!(dst, 1, src, 1)
end

function Base.copyto!(dst::MutableBioSequence{A}, doff::Integer,
                    src::MutableBioSequence{A}, soff::Integer) where {A}
    return copyto!(dst, doff, src, soff, length(src) - soff + 1)
end

function Base.copyto!(dst::MutableBioSequence{A}, doff::Integer,
                    src::MutableBioSequence{A}, soff::Integer, len::Integer) where {A}
    checkbounds(dst, doff:doff+len-1)
    checkbounds(src, soff:soff+len-1)

    if dst.shared || (dst === src && doff > soff)
        orphan!(dst, length(dst), true)
    end

    id = BitIndex(dst, doff)
    is = BitIndex(src, soff)
    rest = len * bitsof(A)

    while rest > 0
        # move `k` bits from `src` to `dst`
        x = dst.data[index(id)]
        y = src.data[index(is)]
        if offset(id) < offset(is)
            y >>= offset(is) - offset(id)
            k = min(64 - offset(is), rest)
        else
            y <<= offset(id) - offset(is)
            k = min(64 - offset(id), rest)
        end
        m = bitmask(k) << offset(id)
        dst.data[index(id)] = y & m | x & ~m

        id += k
        is += k
        rest -= k
    end

    return dst
end

# Actually, users don't need to create a copy of a sequence.
function Base.copy(seq::MutableBioSequence{A}) where {A}
    newseq = MutableBioSequence{A}(seq, 1:lastindex(seq))
    orphan!(newseq, length(seq), true)  # force orphan!
    @assert newseq.data !== seq.data
    return newseq
end


# Encoded copying
# ---------------

function encode_copy!(dst::MutableBioSequence{A},
                      src::Union{AbstractVector,AbstractString}) where {A}
    return encode_copy!(dst, 1, src, 1)
end

function encode_copy!(dst::MutableBioSequence{A},
                      doff::Integer,
                      src::Union{AbstractVector,AbstractString},
                      soff::Integer) where {A}
    return encode_copy!(dst, doff, src, soff, length(src) - soff + 1)
end

function encode_copy!(dst::MutableBioSequence{A},
                      doff::Integer,
                      src::Union{AbstractVector,AbstractString},
                      soff::Integer,
                      len::Integer) where {A}
    if soff != 1 && isa(src, AbstractString) && !isascii(src)
        throw(ArgumentError("source offset ≠ 1 is not supported for non-ASCII string"))
    end

    checkbounds(dst, doff:doff+len-1)
    if length(src) < soff + len - 1
        throw(ArgumentError("source string does not contain $len elements from $soff"))
    end

    orphan!(dst)
    next = BitIndex(dst, doff)
    stop = BitIndex(dst, doff + len)
    i = soff
    while next < stop
        x = UInt64(0)
        j = index(next)
        while index(next) == j && next < stop
            char = src[i]
            i = nextind(src, i)
            x |= enc64(dst, convert(Char, char)) << offset(next)
            next += bitsof(A)
        end
        dst.data[j] = x
    end
    return dst
end

function encode_copy!(dst::MutableBioSequence{A}, doff::Integer,
                      src::AbstractVector{UInt8}, soff::Integer, len::Integer) where {A<:Union{DNAAlphabet{4},RNAAlphabet{4}}}
    checkbounds(dst, doff:doff+len-1)
    if length(src) < soff + len - 1
        throw(ArgumentError("source string does not contain $len elements from $soff"))
    end

    orphan!(dst)
    charmap = A <: DNAAlphabet ? BioSymbols.char_to_dna : BioSymbols.char_to_rna
    i = soff
    next = BitIndex(dst, doff)
    stop = BitIndex(dst, doff + len)

    # head
    if offset(next) != 0
        for d in 0:div(64 - offset(next), 4)-1
            dst[doff+d] = charmap[src[i+d]+1]
        end
        i += div(64 - offset(next), 4)
        next += 64 - offset(next)
    end

    # body
    D = 16
    while next < (stop - offset(stop))
        x::UInt64 = 0
        check = 0x00
        @inbounds for d in 0:D-1
            y = reinterpret(UInt8, charmap[src[i+d]+1])
            x |= UInt64(y) << 4d
            check |= y
        end
        if check & 0x80 != 0
            # invalid byte(s) is detected
            for d in 0:D-1
                if !isvalid(charmap[src[i+d]+1])
                    error("cannot encode $(src[i+d])")
                end
            end
        end
        dst.data[index(next)] = x
        i += D
        next += 64
    end

    # tail
    for d in 0:div(stop - next, 4)-1
        dst[doff+i-soff+d] = charmap[src[i+d]+1]
    end

    return dst
end
