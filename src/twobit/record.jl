# 2bit Record
# ===========

mutable struct Record
    filled::Bool
    dnasize::UInt32
    blockcount::UInt32
    blockstarts::Vector{UInt32}
    blocksizes::Vector{UInt32}
    maskedblockcount::UInt32
    maskedblockstarts::Vector{UInt32}
    maskedblocksizes::Vector{UInt32}
    reserved::UInt32
    packeddna::Vector{UInt8}
end

"""
    TwoBit.Record()

Create an unfilled 2bit record.
"""
function Record()
    return Record(
        false,
        # dnasize-blocksizes
        0, 0, UInt32[], UInt32[],
        # maskedblockcount-reserved
        0, UInt32[], UInt32[], 0,
        # packeddna
        UInt8[])
end

function initialize!(record::Record)
    record.filled = false
    record.dnasize = 0
    record.blockcount = 0
    empty!(record.blockstarts)
    empty!(record.blocksizes)
    record.maskedblockcount = 0
    empty!(record.maskedblockstarts)
    empty!(record.maskedblocksizes)
    record.reserved = 0
    return record
end

function isfilled(record::Record)
    return record.filled
end

function Base.:(==)(record1::Record, record2::Record)
    if isfilled(record1) == isfilled(record2) == true
        return (
            record1.dnasize           == record2.dnasize           &&
            record1.blockcount        == record2.blockcount        &&
            record1.blockstarts       == record2.blockstarts       &&
            record1.blocksizes        == record2.blocksizes        &&
            record1.maskedblockcount  == record2.maskedblockcount  &&
            record1.maskedblockstarts == record2.maskedblockstarts &&
            record1.maskedblocksizes  == record2.maskedblocksizes  &&
            record1.reserved          == record2.reserved          &&
            memcmp(pointer(record1.packeddna), pointer(record2.packeddna), cld(record1.dnasize, 4)) == 0)
    else
        return isfilled(record1) == isfilled(record2) == false
    end
end

function Base.show(io::IO, record::Record)
    print(io, summary(record), ':')
    if isfilled(record)
        println(io)
        println(io, "         length: ", record.dnasize)
        println(io, "       sequence: ", truncate(sequence(record), 40))
        println(io, "       N blocks: ", record.blockcount)
          print(io, "  masked blocks: ", record.maskedblockcount)
    else
        print(io, " <not filled>")
    end
end

function truncate(seq, len)
    if length(seq) > len
        return "$(seq[1:len - 1])…"
    else
        return string(seq)
    end
end


# Accessors
# ---------

function hassequence(record::Record)
    return isfilled(record)
end

"""
    sequence([::Type{S},] record::Record)::S

Get the sequence of `record` as `S`.

`S` can be either `BioSequences.ReferenceSequence` or `BioSequences.DNASequence`.
If `S` is omitted, the default type is `BioSequences.ReferenceSequence`.
"""
function sequence(record::Record)
    return sequence(BioSequences.ReferenceSequence, record)
end

function sequence(::Type{BioSequences.ReferenceSequence}, record::Record)
    checkfilled(record)
    data = decode_sequence(record.packeddna, record.dnasize, BioSequences.BitsPerSymbol{2}(), twobit2refseq_table)
    nmask = falses(record.dnasize)
    for i in 1:record.blockcount
        nmask[record.blockstarts[i] .+ (1:record.blocksizes[i])] .= true
    end
    return BioSequences.ReferenceSequence(data, BioSequences.NMask(nmask), 1:record.dnasize)
end

function sequence(::Type{BioSequences.DNASequence}, record::Record)
    checkfilled(record)
    data = decode_sequence(record.packeddna, record.dnasize, BioSequences.BitsPerSymbol{4}(), twobit2dnaseq_table)
    seq = BioSequences.DNASequence(data, 1:record.dnasize, false)
    for i in 1:record.blockcount
        for j in record.blockstarts[i] .+ (1:record.blocksizes[i])
            seq[j] = BioSequences.DNA_N
        end
    end
    return seq
end

function BioCore.sequence(record::Record)
    return sequence(record)
end

function BioCore.hassequence(record::Record)
    return hassequence(record)
end

"""
    maskedblocks(record::Record)::Vector{UnitRange{Int}}

Get the masked blocks.
"""
function maskedblocks(record::Record)
    checkfilled(record)
    blocks = Vector{UnitRange{Int}}(undef, record.maskedblockcount)
    for i in 1:lastindex(blocks)
        blocks[i] = record.maskedblockstarts[i] .+ (1:record.maskedblocksizes[i])
    end
    return blocks
end

function decode_sequence(packeddna::Vector{UInt8}, seqlen::UInt32, nbits::BioSequences.BitsPerSymbol{n}, table::Vector{UInt64}) where {n}
    @assert n ∈ (2, 4)
    data = zeros(UInt64, cld(seqlen, div(64, n)))
    stop = BioSequences.bitindex(nbits, UInt64, seqlen)
    i = BioSequences.bitindex(nbits, UInt64, 1)
    j = 1
    while i ≤ stop
        data[BioSequences.index(i)] |= table[Int(packeddna[j]) + 1] << BioSequences.offset(i)
        i += 4 * n
        j += 1
    end
    return data
end

const twobit2refseq_table = let
    # T: 00, C: 01, A: 10, G: 11
    f(x) = x == 0b00 ? UInt64(3) :
           x == 0b01 ? UInt64(1) :
           x == 0b10 ? UInt64(0) :
           x == 0b11 ? UInt64(2) : error()
    tcag = 0b00:0b11
    tbl = UInt64[]
    for x in tcag, y in tcag, z in tcag, w in tcag
        push!(tbl, f(x) | f(y) << 2 | f(z) << 4 | f(w) << 6)
    end
    tbl
end

const twobit2dnaseq_table = let
    # T: 00, C: 01, A: 10, G: 11
    f(x) = x == 0b00 ? convert(UInt64, BioSequences.DNA_T) :
           x == 0b01 ? convert(UInt64, BioSequences.DNA_C) :
           x == 0b10 ? convert(UInt64, BioSequences.DNA_A) :
           x == 0b11 ? convert(UInt64, BioSequences.DNA_G) : error()
    tcag = 0b00:0b11
    tbl = UInt64[]
    for x in tcag, y in tcag, z in tcag, w in tcag
        push!(tbl, f(x) | f(y) << 4 | f(z) << 8 | f(w) << 12)
    end
    tbl
end

function checkfilled(record)
    if !isfilled(record)
        throw(ArgumentError("unfilled 2bit record"))
    end
end

function memcmp(p1::Ptr, p2::Ptr, n::Integer)
    return ccall(:memcmp, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Csize_t), p1, p2, n)
end

"""
`WriteRecord{S,T}` is a type holding a named sequence of type `S`, along with
a mask.

Mostly to make writing 2bit sequences painless with the depreceation of
SeqRecord.
I can't see this existing much beyond the addition of masked sequences.
Creates a sequence record suitable for writing to 2bit, in a similar way that
FASTA and FASTQ records are created before being written to file. i.e. by calling
a `Record` method on a name, some sequence, and some masks.
"""
mutable struct WriteRecord{S<:BioSequences.BioSequence}
    name::String
    seq::S
    masks::Union{Vector{UnitRange{Int}}, Nothing}
end

"""
Record(name::AbstractString, seq::BioSequences.BioSequence, masks::Union{Vector{UnitRange{Int}}, Nothing} = nothing)

Prepare a record for writing to a 2bit formatted file.

Needs a `name`, a `sequence`, and (optionally) `masks`: a vector of
ranges that delineate masked regions of sequence.
"""
function Record(name::AbstractString,
                seq::BioSequences.BioSequence,
                masks::Union{Vector{UnitRange{Int}}, Nothing} = nothing)
    return WriteRecord(string(name), seq, masks)
end
