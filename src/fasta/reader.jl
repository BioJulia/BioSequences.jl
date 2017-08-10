# FASTA Reader
# ============

struct Reader{T<:Union{String,IO}} <: BioCore.IO.AbstractReader
    source::T
    index::Nullable{Index}

    # This field is for backward compatibility.
    stream::TranscodingStreams.TranscodingStream

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
    return Reader{typeof(input)}(input, index, TranscodingStreams.NoopStream(input))
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
    filepath = convert(String, input)
    if index == :fai
        faipath = string(filepath, ".fai")
        index = isfile(faipath) ? faipath : nothing
    end
    if index isa Index || index isa Nullable{Index} || index == nothing
        # ok
    elseif index isa AbstractString
        index = Index(index)
    else
        throw(ArgumentError("invalid index argument"))
    end
    return Reader{String}(filepath, index)
end

function Base.eltype(::Type{<:Reader})
    return Record
end

function BioCore.IO.eachrecord(reader::Reader{String}; copy::Bool=true)
    if endswith(reader.source, ".gz")
        stream = CodecZlib.GzipDecompressionStream(open(reader.source))
    else
        stream = open(reader.source)
    end
    return Iterator(stream, copy, true)
end

function BioCore.IO.eachrecord(reader::Reader{<:IO}; copy::Bool=true)
    return Iterator(stream.source, copy, false)
end

function Base.getindex(reader::Reader{String}, name::AbstractString)
    return Reader(open(reader), index=reader.index)[name]
end

function Base.getindex(reader::Reader{<:IO}, name::AbstractString)
    if isnull(reader.index)
        throw(ArgumentError("no index attached"))
    end
    seekrecord(reader.source, get(reader.index), name)
    return first(reader)
end

function Base.start(reader::Reader{String})
    iter = Iterator(open(reader.source), #=copy=#true, #=close=#true)
    state = start(iter)
    return iter, state
end

function Base.start(reader::Reader{<:IO})
    iter = Iterator(reader.source, #=copy=#true, #=close=#false)
    state = start(iter)
    return iter, state
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
    state = IteratorState(1, false, false)
    readrecord!(reader.stream, record, state)
    if !state.read
        throw(ArgumentError("failed to read a FASTA record"))
    end
    return record
end


# Iterator
# --------

struct Iterator
    # input stream
    stream::TranscodingStreams.TranscodingStream

    # return a copy?
    copy::Bool

    # close stream at the end?
    close::Bool

    # placeholder
    record::Record

    function Iterator(stream::TranscodingStreams.TranscodingStream, copy::Bool, close::Bool)
        return new(stream, copy, close, Record())
    end
end

function Iterator(stream::IO, copy::Bool, close::Bool)
    # Wrap an I/O stream with NoopStream.
    return Iterator(TranscodingStreams.NoopStream(stream), copy, close)
end

function Base.iteratorsize(::Type{Iterator})
    return Base.SizeUnknown()
end

function Base.eltype(::Type{Iterator})
    return Record
end

mutable struct IteratorState
    # the current line number
    linenum::Int

    # read a new record?
    read::Bool

    # consumed all input?
    done::Bool
end

function Base.start(iter::Iterator)
    return IteratorState(1, false, false)
end

function Base.done(iter::Iterator, state::IteratorState)
    if state.done
        return true
    elseif !state.read
        readrecord!(iter.stream, iter.record, state)
        if iter.close && state.done
            close(iter.stream)
        end
    end
    return !state.read
end

function Base.next(iter::Iterator, state::IteratorState)
    @assert state.read
    record = iter.record
    if iter.copy
        record = copy(record)
    end
    state.read = false
    return record, state
end

function index!(record::Record)
    stream = TranscodingStreams.NoopStream(IOBuffer(record.data))
    state = IteratorState(1, false, false)
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
            Record, machine, actions, initcode, exitcode, loopunroll=0))
end
