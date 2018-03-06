# FASTQ Reader
# ============

struct Reader{T<:Union{AbstractString,IO}} <: BioCore.IO.AbstractReader
    # data source
    source::T

    # sequencee transformer
    transform::Union{Function,Void}
end

"""
    FASTQ.Reader(input::IO; fill_ambiguous=nothing)

Create a data reader of the FASTQ file format.

# Arguments
* `input`: data source
* `fill_ambiguous=nothing`: fill ambiguous symbols with the given symbol
"""
function Reader(input::Union{AbstractString,IO}; fill_ambiguous=nothing)
    if fill_ambiguous === nothing
        transform = nothing
    else
        transform = generate_fill_ambiguous(fill_ambiguous)
    end
    return Reader{typeof(input)}(input, transform)
end

function Base.eltype(::Type{<:Reader})
    return Record
end

function BioCore.IO.stream(reader::Reader)
    return reader.source
end

function eachrecord(reader::Reader{<:IO})
    return RecordIterator(reader.source, reader.transform)
end

function eachrecord(reader::Reader{<:AbstractString})
    return RecordIterator(open(reader.source), reader.transform)
end

function Base.start(reader::Reader)
    iter = eachrecord(reader)
    return (iter, start(iter))
end

function Base.done(::Reader, state)
    return done(state[1], state[2])
end

function Base.next(::Reader, state)
    return next(state[1], state[2])[1], state
end

#function Base.read(reader::Reader{<:IO}, record::Record)
#end

function Base.close(reader::Reader)
    if reader.source isa IO
        close(reader.source)
    end
    return nothing
end

struct RecordIterator
    stream::TranscodingStream
    transform::Union{Function,Void}

    function RecordIterator(stream::IO, transform)
        if !(stream isa TranscodingStream)
            stream = TranscodingStreams.NoopStream(stream)
        end
        return new(stream, transform)
    end
end

mutable struct RecordIteratorState
    # machine state
    state::Int
    # line number
    linenum::Int
    # is record filled?
    filled::Bool
    # placeholder
    record::Record
end

function Base.iteratorsize(::Type{RecordIterator})
    return Base.SizeUnknown()
end

function Base.eltype(::Type{RecordIterator})
    return Record
end

function Base.start(iter::RecordIterator)
    return RecordIteratorState(1, 1, false, Record())
end

function Base.done(iter::RecordIterator, state::RecordIteratorState)
    if state.filled
        return true
    end
    state.state, state.linenum, state.filled = readrecord!(iter.stream, state.record, (state.state, state.linenum), iter.transform)
    if !state.filled && state.state != 0
        throw(ArgumentError("malformed FASTQ file"))
    end
    return !state.filled
end

function Base.next(iter::RecordIterator, state::RecordIteratorState)
    @assert state.filled
    record = copy(state.record)
    state.filled = false
    return record, state
end

function generate_fill_ambiguous(symbol::BioSymbols.DNA)
    certain = map(UInt8, ('A', 'C', 'G', 'T', 'a', 'c', 'g', 't'))
    # return transform function
    return function (data, range)
        fill = convert(UInt8, convert(Char, symbol))
        for i in range
            if data[i] ∉ certain
                data[i] = fill
            end
        end
        return data
    end
end

machine = (function ()
    re = Automa.RegExp

    hspace = re"[ \t\v]"

    header1 = let
        identifier = re.rep(re.any() \ re.space())
        identifier.actions[:enter] = [:pos]
        identifier.actions[:exit]  = [:header1_identifier]

        description = re.cat(re.any() \ hspace, re"[^\r\n]*")
        description.actions[:enter] = [:pos]
        description.actions[:exit]  = [:header1_description]

        re.cat('@', identifier, re.opt(re.cat(re.rep1(hspace), description)))
    end

    sequence = re"[A-z]*"
    sequence.actions[:enter] = [:pos]
    sequence.actions[:exit]  = [:sequence]

    header2 = let
        identifier = re.rep1(re.any() \ re.space())

        description = re.cat(re.any() \ hspace, re"[^\r\n]*")

        re.cat('+', re.opt(re.cat(identifier, re.opt(re.cat(re.rep1(hspace), description)))))
    end
    header2.actions[:enter] = [:pos]
    header2.actions[:exit]  = [:header2]

    quality = re"[!-~]*"
    quality.actions[:enter] = [:pos]
    quality.actions[:exit]  = [:quality]

    newline = let
        lf = re"\n"
        lf.actions[:enter] = [:countline]

        re.cat(re.opt('\r'), lf)
    end

    record = re.cat(header1, newline, sequence, newline, header2, newline, quality)
    record.actions[:enter] = [:mark]
    record.actions[:exit] = [:record]

    fastq = re.rep(re.cat(record, newline))

    Automa.compile(fastq)
end)()

#write("fastq.dot", Automa.machine2dot(machine))
#run(`dot -Tsvg -o fastq.svg fastq.dot`)

function appendfrom!(dst, dpos, src, spos, n)
    if length(dst) < dpos + n - 1
        resize!(dst, dpos + n - 1)
    end
    copy!(dst, dpos, src, spos, n)
    return dst
end

function is_same_mem(data, pos1, pos2, len)
    checkbounds(data, 1:pos1+len-1)
    checkbounds(data, 1:pos2+len-1)
    return ccall(:memcmp, Cint, (Ptr{Void}, Ptr{Void}, Csize_t), pointer(data, pos1), pointer(data, pos2), len) == 0
end

actions = Dict(
    :mark => :(@mark),
    :pos => :(pos = @relpos(p)),
    :countline => :(linenum += 1),
    :header1_identifier => :(record.identifier = pos:@relpos(p-1)),
    :header1_description => :(record.description = pos:@relpos(p-1)),
    :sequence => :(record.sequence = pos:@relpos(p-1)),
    :header2 => :(second_header_pos = pos+1; second_header_len = @relpos(p)-(pos+1)),
    :quality => :(record.quality = pos:@relpos(p-1)),
    :record => quote
        appendfrom!(record.data, 1, data, @markpos, p-@markpos)
        record.filled = 1:(p-@markpos)
        found = true
        @escape
    end,
)
initcode = quote
    pos = 0
    second_header_pos = 0
    second_header_len = 0
    found = false
    initialize!(record)
    cs, linenum = state
end
loopcode = quote
    if cs < 0
        throw(ArgumentError("malformed FASTQ file at line $(linenum)"))
    elseif found && length(record.sequence) != length(record.quality)
        throw(ArgumentError("mismatched sequence and quality length"))
    elseif found && second_header_len > 0
        first_header_pos = first(record.identifier)
        first_header_len = max(last(record.identifier), last(record.description)) - first_header_pos + 1
        if first_header_len != second_header_len || !is_same_mem(record.data, first_header_pos, second_header_pos, first_header_len)
            throw(ArgumentError("mismatched headers"))
        end
    elseif found && transform != nothing
        transform(record.data, record.sequence)
    end
    found && @goto __return__
end
returncode = :(return cs, linenum, found)
context = Automa.CodeGenContext(generator=:goto)
Automa.Stream.generate_reader(
    :readrecord!,
    machine,
    arguments=(:(record::Record),:(state::Tuple{Int,Int}),:(transform)),
    actions=actions,
    context=context,
    initcode=initcode,
    loopcode=loopcode,
    returncode=returncode,
) |> eval

function index!(record::Record)
    stream = TranscodingStreams.NoopStream(IOBuffer(record.data))
    cs, linenum, found = readrecord!(stream, record, (1, 1), nothing)
    if !found || !allspace(stream)
        throw(ArgumentError("invalid FASTQ record"))
    end
    return record
end

function allspace(stream)
    while !eof(stream)
        if !isspace(read(stream, Char))
            return false
        end
    end
    return true
end
