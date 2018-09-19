# ReferenceSequence
# =================
#
# DNA sequence for reference genomes.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

"""
Reference Sequence

Reference sequence is a sequence of A/C/G/T/N. In the internals, it compresses
`N` positions and consumes less than three bits per base. Unlike `BioSequence`,
reference sequences are immutable and hence no modifyting operators are
provided.
"""
struct ReferenceSequence <: BioSequence
    data::Vector{UInt64}  # 2-bit encoding of A/C/G/T nucloetides
    nmask::NMask          # positions of N
    part::UnitRange{Int}  # interval within `data` defining the (sub)sequence
end

Base.length(seq::ReferenceSequence) = length(seq.part)
Alphabet(::Type{ReferenceSequence}) = DNAAlphabet{2}()
Base.summary(seq::ReferenceSequence) = string(length(seq), "nt Reference Sequence")

function Base.copy(seq::ReferenceSequence)
    return ReferenceSequence(copy(seq.data), copy(seq.nmask), seq.part)
end

function ReferenceSequence()
    return ReferenceSequence(UInt64[], NMask(), 1:0)
end

function ReferenceSequence(seq::ReferenceSequence, part::UnitRange{<:Integer})
    ReferenceSequence(seq.data, seq.nmask, part)
end

function ReferenceSequence(src::Vector{UInt8}, startpos::Integer=1,
                           len::Integer=length(src))
    return encode(src, startpos, len)
end

function ReferenceSequence(seq::GeneralSequence{<:DNAAlphabet})
    data = Vector{UInt64}(undef, cld(length(seq), 32))
    nmask = falses(length(seq))
    i = 1
    for j in 1:lastindex(data)
        x = UInt64(0)
        r = 0
        while r < 64 && i ≤ lastindex(seq)
            nt = seq[i]
            if nt == DNA_A
                x |= convert(UInt64, 0) << r
            elseif nt == DNA_C
                x |= convert(UInt64, 1) << r
            elseif nt == DNA_G
                x |= convert(UInt64, 2) << r
            elseif nt == DNA_T
                x |= convert(UInt64, 3) << r
            elseif nt == DNA_N
                nmask[i] = true
            else
                throw(ArgumentError("invalid symbol $(seq[i]) ∉ {A,C,G,T,N} at $i"))
            end
            i += 1
            r += 2
        end
        data[j] = x
    end
    return ReferenceSequence(data, NMask(nmask), 1:length(seq))
end

function ReferenceSequence(str::AbstractString)
    if !isascii(str)
        throw(ArgumentError("attempt to convert a non-ASCII string to ReferenceSequence"))
    end
    return encode([UInt8(char) for char in str], 1, length(str))
end

function DNASequence(seq::ReferenceSequence)
    bioseq = DNASequence(length(seq))
    for i in 1:lastindex(seq)
        bioseq[i] = seq[i]
    end
    return bioseq
end

function Base.convert(::Type{S}, seq::ReferenceSequence) where {S<:AbstractString}
    return S([Char(nt) for nt in seq])
end
Base.String(seq::ReferenceSequence) = convert(String, seq)

@inline function bitindex(seq::ReferenceSequence, i::Integer)
    return bitindex(BitsPerSymbol{2}(), UInt64, i + first(seq.part) - 1)
end

# Create ReferenceSequence object from the ascii-encoded `data`
function encode(src::Vector{UInt8}, from::Integer, len::Integer)
    #println("src: ", src)
    data = zeros(UInt64, cld(len, 32))
    nmask = falses(len)
    #next = bitindex(1, 2)
    #next = BitIndex{2,UInt64}(1)
    next = bitindex(BitsPerSymbol{2}(), UInt64, 1)
    #stop = bitindex(len + 1, 2)
    #stop = BitIndex{2, UInt64}(len + 1)
    stop = bitindex(BitsPerSymbol{2}(), UInt64, len + 1)
    #println("next: ", next)
    #println("stop: ", stop)
    i = from
    while next < stop
        x = UInt64(0)
        j = index(next)
        #println("x: ", x)
        #println("j: ", j)
        while index(next) == j && next < stop
            # FIXME: Hotspot
            char = convert(Char, src[i])
            nt = convert(DNA, char)
            #println("char: ", char)
            #println("nt: ", nt)
            #println("!isambiguous: ", !isambiguous(nt))
            if !isambiguous(nt)
                #println("Encoded nt: ", hex(encode(DNAAlphabet{2}, nt) << offset(next)))
                x |= UInt64(encode(DNAAlphabet{2}(), nt)) << offset(next)
            elseif nt == DNA_N
                nmask[i] = true
            else
                error("'", char, "'", " is not allowed")
            end
            i += 1
            next += 2
            #println("i and next: ", i, ", ", next)
        end
        data[j] = x
    end
    return ReferenceSequence(data, NMask(nmask), 1:len)
end

function Base.checkbounds(seq::ReferenceSequence, part::UnitRange)
    if isempty(part) || (1 ≤ first(part) && last(part) ≤ lastindex(seq))
        return true
    end
    throw(BoundsError(seq, part))
end

@inline function inbounds_getindex(seq::ReferenceSequence, i::Integer)
    if seq.nmask[i + first(seq.part) - 1]
        return DNA_N
    else
        j = bitindex(seq, i)
        return DNA(0x01 << ((seq.data[index(j)] >> offset(j)) & 0b11))
    end
end

function Base.getindex(seq::ReferenceSequence, part::UnitRange{<:Integer})
    checkbounds(seq, part)
    return ReferenceSequence(seq, part)
end

function find_next_ambiguous(seq::ReferenceSequence, i::Integer)
    return findnextn(seq.nmask, i)
end
