# 2bit Writer
# ===========

mutable struct Writer{T<:IO} <: BioCore.IO.AbstractWriter
    # output stream
    output::T
    # sequence names
    names::Vector{String}
    # bit vector to check if each sequence is already written or not
    written::BitVector
end

"""
    TwoBitWriter(output::IO, names::AbstractVector)

Create a data writer of the 2bit file format.

# Arguments
* `output`: data sink
* `names`: a vector of sequence names written to `output`
"""
function Writer(output::IO, names::AbstractVector)
    writer = Writer(output, names, falses(length(names)))
    write_header(writer)
    write_index(writer)
    return writer
end

function BioCore.IO.stream(writer::Writer)
    return writer.output
end

function Base.close(writer::Writer)
    if !all(writer.written)
        error("one or more sequences are not written")
    end
    close(writer.output)
end

function write_header(writer::Writer)
    n = 0
    n += write(writer.output, SIGNATURE)
    n += write(writer.output, UInt32(0))
    n += write(writer.output, UInt32(length(writer.names)))
    n += write(writer.output, UInt32(0))
    return n
end

function write_index(writer::Writer)
    n = 0
    for name in writer.names
        n += write(writer.output, UInt8(length(name)))
        n += write(writer.output, name)
        n += write(writer.output, UInt32(0))  # filled later
    end
    return n
end

# Update the file offset of a sequence in the index section.
function update_offset(writer::Writer, seqname, seqoffset)
    @assert seqname ∈ writer.names
    offset = 16
    for name in writer.names
        offset += sizeof(UInt8) + length(name)
        if name == seqname
            old = position(writer.output)
            seek(writer.output, offset)
            write(writer.output, UInt32(seqoffset))
            seek(writer.output, old)
            return
        end
        offset += sizeof(UInt32)
    end
end

function Base.write(writer::Writer, record::WriteRecord)
    i = findfirst(isequal(record.name), writer.names)
    if i == 0
        error("sequence \"", record.name, "\" doesn't exist in the writing list")
    elseif writer.written[i]
        error("sequence \"", record.name, "\" is already written")
    end

    output = writer.output
    update_offset(writer, record.name, position(output))

    n = 0
    n += write(output, UInt32(length(record.seq)))
    n += write_n_blocks(output, record.seq)
    n += write_masked_blocks(output, record.masks)
    n += write(output, UInt32(0))  # reserved bytes
    n += write_twobit_sequence(output, record.seq)

    writer.written[i] = true
    return n
end

function make_n_blocks(seq)
    starts = UInt32[]
    sizes = UInt32[]
    i = 1
    while i ≤ lastindex(seq)
        nt = seq[i]
        if nt == BioSequences.DNA_N
            start = i - 1  # 0-based index
            push!(starts, start)
            while i ≤ lastindex(seq) && seq[i] == BioSequences.DNA_N
                i += 1
            end
            push!(sizes, (i - 1) - start)
        elseif BioSequences.isambiguous(nt)
            error("ambiguous nucleotide except N is not supported in the 2bit file format")
        else
            i += 1
        end
    end
    return starts, sizes
end

function write_n_blocks(output, seq)
    blockstarts, blocksizes = make_n_blocks(seq)
    @assert length(blockstarts) == length(blocksizes)
    n = 0
    n += write(output, UInt32(length(blockstarts)))
    n += write(output, blockstarts)
    n += write(output, blocksizes)
    return n
end

function write_masked_blocks(output, masks)
    n = 0
    if masks != nothing
        n += write(output, UInt32(length(masks)))
        for mblock in masks
            n += write(output, UInt32(first(mblock) - 1))  # 0-based
        end
        for mblock in masks
            n += write(output, UInt32(length(mblock)))
        end
    else
        n += write(output, UInt32(0))
    end
    return n
end

function write_twobit_sequence(output, seq)
    n = 0
    i = 4
    while i ≤ lastindex(seq)
        x::UInt8 = 0
        x |= nuc2twobit(seq[i - 3]) << 6
        x |= nuc2twobit(seq[i - 2]) << 4
        x |= nuc2twobit(seq[i - 1]) << 2
        x |= nuc2twobit(seq[i - 0]) << 0
        n += write(output, x)
        i += 4
    end
    r = length(seq) % 4
    if r > 0
        let x::UInt8 = 0
            i = lastindex(seq) - r + 1
            while i ≤ lastindex(seq)
                x = x << 2 | nuc2twobit(seq[i])
                i += 1
            end
            x <<= (4 - r) * 2
            n += write(output, x)
        end
    end
    return n
end

function nuc2twobit(nt::BioSequences.DNA)
    return (
        nt == BioSequences.DNA_A ? 0b10 :
        nt == BioSequences.DNA_C ? 0b01 :
        nt == BioSequences.DNA_G ? 0b11 :
        nt == BioSequences.DNA_T ? 0b00 :
        nt == BioSequences.DNA_N ? 0b00 : error()
    )
end
