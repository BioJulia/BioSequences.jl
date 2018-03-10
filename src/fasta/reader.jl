# FASTA Reader
# ============

struct Reader{T<:Union{IO,AbstractString}} <: BioCore.IO.AbstractReader
    source::T
    index::Nullable{Index}
end

"""
    FASTA.Reader(input::IO; index=nothing)

Create a data reader of the FASTA file format.

# Arguments
* `input`: data source
* `index=nothing`: filepath to a random access index (currently *fai* is supported)
"""
function Reader(input::Union{AbstractString,IO}; index=nothing)
    if isa(index, AbstractString)
        index = Index(index)
    else
        if index != nothing
            throw(ArgumentError("index must be a filepath or nothing"))
        end
    end
    return Reader{typeof(input)}(input, index)
end

function Base.eltype(::Type{<:Reader})
    return Record
end

function BioCore.IO.stream(reader::Reader)
    return reader.source
end

function ReaderIterTools.eachrecord(reader::Reader{<:IO})
    return RecordIterator(reader.source)
end

function ReaderIterTools.eachrecord(reader::Reader{<:AbstractString})
    return RecordIterator(open(reader.source))
end

function Base.start(reader::Reader)
    iter = ReaderIterTools.eachrecord(reader)
    return (iter, start(iter))
end

function Base.done(::Reader, state)
    return done(state[1], state[2])
end

function Base.next(::Reader, state)
    return next(state[1], state[2])[1], state
end

function Base.close(reader::Reader)
    if reader.source isa IO
        close(reader.source)
    end
    return nothing
end

function Base.getindex(reader::Reader, name::AbstractString)
    if isnull(reader.index)
        throw(ArgumentError("no index attached"))
    end
    seekrecord(reader.source, get(reader.index), name)
    record = Record()
    cs, linenum, found = readrecord!(TranscodingStreams.NoopStream(reader.source), record, (1, 1))
    @assert cs ≥ 0 && found
    return record
end

struct RecordIterator <: ReaderIterTools.AbstractRecordIterator{Record}
    stream::TranscodingStream

    function RecordIterator(stream::IO)
        if !(stream isa TranscodingStream)
            stream = TranscodingStreams.NoopStream(stream)
        end
        return new(stream)
    end
end

function Base.done(iter::RecordIterator, state)
    if state.filled
        return true
    end
    state.state, state.linenum, state.filled = readrecord!(iter.stream, state.record, (state.state, state.linenum))
    if !state.filled && state.state != 0
        throw(ArgumentError("malformed FASTQ file"))
    end
    return !state.filled
end

machine = (function ()
    re = Automa.RegExp

    hspace = re"[ \t\v]"

    identifier = re.rep(re.any() \ re.space())
    identifier.actions[:enter] = [:pos]
    identifier.actions[:exit]  = [:identifier]

    description = re.cat(re.any() \ hspace, re"[^\r\n]*")
    description.actions[:enter] = [:pos]
    description.actions[:exit]  = [:description]

    header = re.cat('>', identifier, re.opt(re.cat(re.rep1(hspace), description)))
    header.actions[:enter] = [:mark]
    header.actions[:exit]  = [:header]

    # '*': terminal, `-': gap
    letters = re"[A-Za-z*\-]+"
    letters.actions[:enter] = [:mark]
    letters.actions[:exit]  = [:letters]

    newline = let
        lf = re"\n"
        lf.actions[:enter] = [:countline]

        re.cat(re.opt('\r'), lf)
    end

    sequence = re.opt(re.cat(letters, re.rep(re.cat(re.rep1(re.alt(hspace, newline)), letters))))
    #sequence.actions[:enter] = [:sequence_start]

    record = re.cat(header, re.rep1(newline), sequence, re.rep1(newline))
    record.actions[:exit] = [:record]

    record_trailing = re.cat(header, re.rep1(newline), sequence)
    record_trailing.actions[:exit] = [:record]

    fasta = re.cat(re.rep(newline), re.rep(record), re.opt(record_trailing))

    Automa.compile(fasta)
end)()

function appendfrom!(dst, dpos, src, spos, n)
    if length(dst) < dpos + n - 1
        resize!(dst, dpos + n - 1)
    end
    copy!(dst, dpos, src, spos, n)
    return dst
end

actions = Dict(
    :mark => :(@mark),
    :pos => :(pos = @relpos(p)),
    :countline => :(linenum += 1),
    :identifier => :(record.identifier = pos:@relpos(p-1)),
    :description => :(record.description = pos:@relpos(p-1)),
    :header => quote
        let n = p - @markpos
            appendfrom!(record.data, 1, data, @markpos, n)
            filled += n
            appendfrom!(record.data, filled + 1, b"\n", 1, 1)
            filled += 1
        end
    end,
    :letters => quote
        let markpos = @markpos(), n = @relpos(p-1) - @relpos(markpos) + 1
            appendfrom!(record.data, filled + 1, data, markpos, n)
            if isempty(record.sequence)
                record.sequence = filled+1:filled+n
            else
                record.sequence = first(record.sequence):last(record.sequence)+n
            end
            filled += n
        end
    end,
    :record => quote
        record.filled = 1:filled
        found = true
        @escape
    end,
)
initcode = quote
    pos = 0
    filled = 0
    found = false
    initialize!(record)
    cs, linenum = state
end
loopcode = quote
    if cs < 0
        throw(ArgumentError("malformed FASTA file at line $(linenum)"))
    end
    found && @goto __return__
end
returncode = :(return cs, linenum, found)
context = Automa.CodeGenContext(generator=:goto)
Automa.Stream.generate_reader(
    :readrecord!,
    machine,
    arguments=(:(record::Record), :(state::Tuple{Int,Int})),
    actions=actions,
    context=context,
    initcode=initcode,
    loopcode=loopcode,
    returncode=returncode,
) |> eval

function index!(record::Record)
    stream = TranscodingStreams.NoopStream(IOBuffer(record.data))
    cs, linenum, found = readrecord!(stream, record, (1, 1))
    if cs != 0
        throw(ArgumentError("invalid FASTA record"))
    end
    return record
end
