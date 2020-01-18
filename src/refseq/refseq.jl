###
### ReferenceSequence
###
###
### DNA sequence for reference genomes.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

"""
Reference Sequence

Reference sequence is a sequence of A/C/G/T/N. In the internals, it compresses
`N` positions and consumes less than three bits per base. Unlike `BioSequence`,
reference sequences are immutable and hence no modifying operators are
provided.
"""
struct ReferenceSequence <: BioSequence{DNAAlphabet{2}}
    data::Vector{UInt64}  # 2-bit encoding of A/C/G/T nucloetides
    nmask::NMask          # positions of N
    part::UnitRange{Int}  # interval within `data` defining the (sub)sequence
end

Base.length(seq::ReferenceSequence) = last(seq.part) - first(seq.part) + 1

Base.summary(seq::ReferenceSequence) = string(length(seq), "nt Reference Sequence")

function Base.copy(seq::ReferenceSequence)
    return ReferenceSequence(copy(seq.data), copy(seq.nmask), seq.part)
end

function ReferenceSequence()
    return ReferenceSequence(UInt64[], NMask(), 1:0)
end

function ReferenceSequence(seq::ReferenceSequence, part::UnitRange{<:Integer})
    checkbounds(seq, part)
    ReferenceSequence(seq.data, seq.nmask, part)
end

function ReferenceSequence(src::Vector{UInt8}, startpos::Integer=1,
                           len::Integer=length(src))
    return encode(src, startpos, len)
end

function ReferenceSequence(seq::LongSequence{<:DNAAlphabet})
    data = Vector{UInt64}(undef, cld(length(seq), 32))
    nmask = falses(length(seq))
    i = 1
    @inbounds for j in 1:lastindex(data)
        x = UInt64(0)
        r = 0
        while (r < 64) & (i ≤ lastindex(seq))
            nt = seq[i]
            enc = twobitnucs[reinterpret(UInt8, nt) + 1]
            if enc == 0xff
                if nt == DNA_N
                    enc = 0x00
                    nmask[i] = true
                else
                    throw(ArgumentError("invalid symbol $(nt) ∉ {A,C,G,T,N} at $i"))
                end
            end
            x |= (enc % UInt64) << (r & 63)
            i += 1
            r += 2
        end
        data[j] = x
    end
    return ReferenceSequence(data, NMask(nmask), 1:length(seq))
end

function ReferenceSequence(str::AbstractString)
    return encode([UInt8(char) for char in str], 1, length(str))
end

function ReferenceSequence(str::Union{String, SubString{String}})
    len = ncodeunits(str) - firstindex(str) + 1
    v = GC.@preserve str unsafe_wrap(Vector{UInt8}, pointer(str), len)
    return encode(v, 1, len)
end


function LongDNASeq(seq::ReferenceSequence)
    bioseq = LongSequence{DNAAlphabet{4}}(length(seq))
    @inbounds for i in 1:lastindex(seq)
        bioseq[i] = seq[i]
    end
    return bioseq
end

@inline function bitindex(seq::ReferenceSequence, i::Integer)
    return bitindex(BitsPerSymbol{2}(), UInt64, i + first(seq.part) - 1)
end

@noinline function throw_refseq_encode_err(byte::UInt8)
    throw(error("Cannot encode $byte to reference DNA"))
end

function encode(src::Vector{UInt8}, from::Integer, len::Integer)
    len < 0 && throw(ArgumentError("length cannot be negative"))
    checkbounds(src, from+len-1)
    data = Vector{UInt64}(undef, seq_data_len(DNAAlphabet{2}, len))
    mask = falses(len)
    chunkrem = 32
    chunki = 1
    sh = 0
    chunk = zero(UInt)
    @inbounds for i in from:from+len-1
        if iszero(chunkrem)
            data[chunki] = chunk
            chunki += 1
            chunk = zero(UInt)
            chunkrem = 32
            sh = 0
        end
        byte = src[i]
        enc = stringbyte(DNAAlphabet{2}(), byte)
        if enc == 0x80
            if stringbyte(DNAAlphabet{4}(), byte) == 0x0f
                mask[i - from + 1] = true
                enc = 0x00
            else
                throw_refseq_encode_err(byte)
            end
        end
        chunk |= (enc % UInt64) << (sh & 63)
        chunkrem -= 1
        sh += 2
    end
    @inbounds data[chunki] = chunk
    nmask = NMask(mask)
    ReferenceSequence(data, nmask, 1:len)
end

function Base.checkbounds(seq::ReferenceSequence, part::UnitRange)
    if isempty(part) | ((1 ≤ first(part)) & (last(part) ≤ lastindex(seq)))
        return true
    end
    throw(BoundsError(seq, part))
end

@inline function inbounds_getindex(seq::ReferenceSequence, i::Integer)
    @inbounds if seq.nmask[i + first(seq.part) - 1]
        return DNA_N
    else
        j = bitindex(seq, i)
        chunk = seq.data[index(j)]
        return reinterpret(DNA, 0x01 << ((chunk >> offset(j)) & 0b11))
    end
end

function Base.getindex(seq::ReferenceSequence, part::UnitRange{<:Integer})
    return ReferenceSequence(seq, part)
end

function find_next_ambiguous(seq::ReferenceSequence, i::Integer)
    return findnextn(seq.nmask, i)
end
