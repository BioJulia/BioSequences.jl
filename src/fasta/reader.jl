# FASTA Reader
# ============

# When T = AbstractString, the source field is a filename.
# When T = IO, the source field is an input I/O stream.
struct Reader{T<:Union{AbstractString,IO}} <: BioCore.IO.AbstractReader
    source::T
    index::Nullable{Index}

    # This field is for backward compatibility.
    stream::TranscodingStream

    function Reader{T}(source::T, index) where T
        return new(source, index)
    end

    function Reader{T}(source::T, index, stream) where T
        return new(source, index, stream)
    end
end

"""
    FASTA.Reader(input::IO; index=nothing)

Create a data reader of the FASTA file format.

# Arguments
* `input`: data source
* `index=nothing`: filepath to a random access index (currently *fai* is supported)
"""
function Reader(input::IO; index=nothing)
    if index isa Index || index isa Nullable{Index} || index == nothing
        # ok
    elseif index isa AbstractString
        index = Index(index)
    else
        throw(ArgumentError("invalid index argument"))
    end
    return Reader{typeof(input)}(input, index, NoopStream(input))
end

"""
    FASTA.Reader(input::AbstractString; index=:fai)

Create a data reader of the FASTA file format.

If `index=:fai` it tries to find the fai index file based on the input filepath
(i.e. `string(filepath, ".fai")`).

# Arguments
- `input`: input filepath
- `index=:fai`: random access index
"""
function Reader(input::AbstractString; index=:fai)
    if index == :fai
        faipath = string(input, ".fai")
        index = isfile(faipath) ? faipath : nothing
    end
    if index isa Index || index isa Nullable{Index} || index == nothing
        # ok
    elseif index isa AbstractString
        index = Index(index)
    else
        throw(ArgumentError("invalid index argument"))
    end
    return Reader{typeof(input)}(input, index)
end

function Base.show(io::IO, reader::Reader)
    print(io, summary(reader), "(<source=$(reader.source),index=$(get(reader.index, "null"))>)")
end

function Base.eltype(::Type{<:Reader})
    return Record
end

function BioCore.IO.eachrecord(reader::Reader{<:AbstractString}; copy::Bool=true)
    file = open(reader.source)
    if endswith(reader.source, ".gz")
        stream = CodecZlib.GzipDecompressorStream(file)
    else
        stream = NoopStream(file)
    end
    return RecordIterator{Record}(stream, copy=copy, close=true)
end

function BioCore.IO.eachrecord(reader::Reader{<:IO}; copy::Bool=true)
    return RecordIterator{Record}(stream.source, copy=copy, close=false)
end

function Base.getindex(reader::Reader{<:AbstractString}, name::AbstractString)
    return open(file -> Reader(file, index=reader.index)[name], reader)
end

function Base.getindex(reader::Reader{<:IO}, name::AbstractString)
    if isnull(reader.index)
        throw(ArgumentError("no index attached"))
    end
    seekrecord(reader.source, get(reader.index), name)
    return first(reader)
end

function Base.start(reader::Reader{<:AbstractString})
    iter = RecordIterator{Record}(open(reader.source), copy=true, close=true)
    return iter, start(iter)
end

function Base.start(reader::Reader{<:IO})
    iter = RecordIterator{Record}(reader.source, copy=true, close=false)
    return iter, start(iter)
end

function Base.done(::Reader, state)
    return done(state[1], state[2])
end

function Base.next(::Reader, state)
    item, _ = next(state[1], state[2])
    return item, state
end

function BioCore.IO.stream(reader::Reader{<:IO})
    return reader.stream
end

function Base.read!(reader::Reader{<:IO}, record::Record)
    state = RecordIteratorState()
    readrecord!(reader.stream, record, state)
    if !state.read
        throw(ArgumentError("failed to read a FASTA record"))
    end
    return record
end


# Format
# ------

function index!(record::Record)
    stream = NoopStream(IOBuffer(record.data))
    state = RecordIteratorState()
    readrecord!(stream, record, state)
    if !state.read || !state.done
        throw(ArgumentError("invalid FASTA record"))
    end
    return record
end

let
    # Define FASTA record machine.
    machine = (function ()
        cat = Automa.RegExp.cat
        rep = Automa.RegExp.rep
        rep1 = Automa.RegExp.rep1
        alt = Automa.RegExp.alt
        opt = Automa.RegExp.opt
        any = Automa.RegExp.any
        space = Automa.RegExp.space

        newline = let
            lf = re"\n"
            lf.actions[:enter] = [:countline]

            cat(opt('\r'), lf)
        end

        hspace = re"[ \t\v]"

        identifier = rep(any() \ space())
        identifier.actions[:enter] = [:pos]
        identifier.actions[:exit]  = [:identifier]

        description = cat(any() \ hspace, re"[^\r\n]*")
        description.actions[:enter] = [:pos]
        description.actions[:exit]  = [:description]

        header = cat('>', identifier, opt(cat(rep1(hspace), description)))
        header.actions[:enter] = [:mark]
        header.actions[:exit]  = [:copy]

        # '*': terminal, `-': gap
        letters = re"[A-Za-z*\-]+"
        letters.actions[:enter] = [:mark]
        letters.actions[:exit]  = [:copy]

        sequence = opt(cat(letters, rep(cat(rep1(alt(hspace, newline)), letters))))
        sequence.actions[:enter] = [:seqstart]

        record = cat(header, rep1(newline), sequence, rep1(newline))
        record.actions[:exit] = [:record]

        # lookahead
        record = cat(record, re">?")

        return Automa.compile(record)
    end)()

    #= debug
    write("fasta.dot", Automa.machine2dot(file_machine))
    run(`dot -Tsvg -o fasta.svg fasta.dot`)
    =#

    # Define init, exit and action code.
    initcode = quote
        initialize!(record)
        nfilled = pos = seqstart = 0
    end
    exitcode = quote
        if cs < 0
            if ismarked(stream)
                unmark(stream)
            end
            throw(ArgumentError("unexpected data at line $(state.linenum)"))
        else
            state.done = cs == 0
            if state.read || state.done
                return
            elseif p > p_eof ≥ 0
                if cs == $(machine.start_state)
                    # empty file
                    return
                elseif seqstart > 0
                    # EOF without newline; this is a hacky workaround and should
                    # be solved by supporting EOF in Automa.jl.
                    let len = p - buffer.markpos
                        if length(record.data) < nfilled + len
                            resize!(record.data, nfilled + len)
                        end
                        copy!(record.data, nfilled + 1, buffer.data, buffer.markpos, len)
                        nfilled += len
                    end
                    @unmark
                    record.sequence = seqstart:nfilled
                    record.filled = 1:nfilled
                    state.read = state.done = true
                    return
                else
                    throw(ArgumentError("unexpected EOF at line $(state.linenum)"))
                end
            end
        end
    end
    actions = Dict(
        :countline => :(state.linenum += 1),
        :pos => :(pos = @pos),
        :identifier => :(record.identifier = pos:(@pos)-1),
        :description => :(record.description = pos:(@pos)-1),
        :mark => :(@mark),
        :copy => quote
            let len = p - buffer.markpos
                if length(record.data) < nfilled + len
                    resize!(record.data, nfilled + len)
                end
                copy!(record.data, nfilled + 1, buffer.data, buffer.markpos, len)
                nfilled += len
            end
            @unmark
            if seqstart == 0
                # Insert a newline byte after the header.
                if length(record.data) == nfilled
                    resize!(record.data, length(record.data) + 1)
                end
                record.data[nfilled+1] = UInt8('\n')
                nfilled += 1
            end
        end,
        :seqstart => :(seqstart = nfilled + 1),
        :record => quote
            record.sequence = seqstart:nfilled
            record.filled = 1:nfilled
            state.read = true
            if p ≤ p_end && data[p] == UInt8('>')
                # decrement the reading position to cancel the lookahead
                p -= 1
            end
            @escape
        end)

    # Generate readrecord! function.
    eval(
        BioCore.ReaderHelper.generate_readrecord_function(
            Record, machine, actions, initcode, exitcode, loopunroll=10))
end
