# FASTQ Reader
# ============

struct Reader{T<:Union{String,IO}} <: BioCore.IO.AbstractReader
    source::T

    # This field is for backward compatibility.
    stream::TranscodingStream

    function Reader{T}(source::T) where T
        return new(source)
    end

    function Reader{T}(source::T, stream) where T
        return new(source, stream)
    end
end

"""
    FASTQ.Reader(input::IO)

Create a data reader of the FASTQ file format.

# Arguments
* `input`: data source
"""
function Reader(input::IO; fill_ambiguous=nothing)
    if fill_ambiguous !== nothing
        warn("fill_ambiguous keyword argument is removed; use fill_ambiguous! instead")
    end
    return Reader{typeof(input)}(input, NoopStream(input))
end

function Reader(input::AbstractString; fill_ambiguous=nothing)
    if fill_ambiguous !== nothing
        warn("fill_ambiguous keyword argument is removed; use fill_ambiguous! instead")
    end
    return Reader{String}(convert(String, input))
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
        #throw(ArgumentError("failed to read a FASTQ record"))
        throw(EOFError())
    end
    return record
end


# Iterator
# --------

# FIXME: Remove redundancy since these are almost identical to those of FASTA.
struct Iterator
    # input stream
    stream::TranscodingStream

    # return a copy?
    copy::Bool

    # close stream at the end?
    close::Bool

    # placeholder
    record::Record

    function Iterator(stream::TranscodingStream, copy::Bool, close::Bool)
        return new(stream, copy, close, Record())
    end
end

function Iterator(stream::IO, copy::Bool, close::Bool)
    return Iterator(NoopStream(stream), copy, close)
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
    stream = NoopStream(IOBuffer(record.data))
    state = IteratorState(1, false, false)
    readrecord!(stream, record, state)
    if !state.read || !state.done
        throw(ArgumentError("invalid FASTQ record"))
    end
    return record
end

function check_identical(data1, range1, data2, range2)
    len = length(range1)
    if len != length(range2) ||
            BioCore.Mem.cmp(pointer(data1, first(range1)), pointer(data2, first(range2)), len) != 0
        throw(ArgumentError("header mismatch"))
    end
end

let
    # Define FASTQ record machine.
    # NOTE: This does not support line-wraps within sequence and quality.
    machine = (function ()
        cat = Automa.RegExp.cat
        rep = Automa.RegExp.rep
        rep1 = Automa.RegExp.rep1
        alt = Automa.RegExp.alt
        opt = Automa.RegExp.opt
        any = Automa.RegExp.any
        space = Automa.RegExp.space

        hspace = re"[ \t\v]"

        header1 = let
            identifier = rep(any() \ space())
            identifier.actions[:enter] = [:pos]
            identifier.actions[:exit]  = [:header1_identifier]

            description = cat(any() \ hspace, re"[^\r\n]*")
            description.actions[:enter] = [:pos]
            description.actions[:exit]  = [:header1_description]

            cat('@', identifier, opt(cat(rep1(hspace), description)))
        end

        sequence = re"[A-z]*"
        sequence.actions[:enter] = [:pos]
        sequence.actions[:exit]  = [:sequence]

        header2 = let
            identifier = rep1(any() \ space())
            identifier.actions[:enter] = [:pos]
            identifier.actions[:exit]  = [:header2_identifier]

            description = cat(any() \ hspace, re"[^\r\n]*")
            description.actions[:enter] = [:pos]
            description.actions[:exit]  = [:header2_description]

            cat('+', opt(cat(identifier, opt(cat(rep1(hspace), description)))))
        end

        quality = re"[!-~]*"
        quality.actions[:enter] = [:pos]
        quality.actions[:exit]  = [:quality]

        newline = let
            lf = re"\n"
            lf.actions[:enter] = [:countline]

            cat(opt('\r'), lf)
        end

        record = cat(
            header1,  newline,
            sequence, newline,
            header2,  newline,
            quality,  newline,
        )
        record.actions[:enter] = [:recstart]
        record.actions[:exit]  = [:record]

        # lookahead
        record = cat(record, re"@?")

        return Automa.compile(record)
    end)()

    #= debug
    write("fastq.dot", Automa.machine2dot(machine))
    run(`dot -Tsvg -o fastq.svg fastq.dot`)
    =#

    # Define init, exit and action code.
    initcode = quote
        initialize!(record)
        pos = recstart = 0
        header2_identifier = header2_description = 1:0
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
            elseif cs == $(machine.start_state)
                # empty file
                return
            elseif p > p_eof ≥ 0
                throw(ArgumentError("unexpected EOF at line $(state.linenum)"))
            end
        end
    end
    actions = Dict(
        :pos       => :(pos = @pos),
        :countline => :(state.linenum += 1),
        :header1_identifier  => :(record.identifier  = pos:(@pos)-1),
        :header1_description => :(record.description = pos:(@pos)-1),
        :header2_identifier  => :(header2_identifier  = pos:(@pos)-1),
        :header2_description => :(header2_description = pos:(@pos)-1),
        :sequence => :(record.sequence = pos:(@pos)-1),
        :quality  => :(record.quality  = pos:(@pos)-1),
        :recstart => :(@mark),
        :record   => quote
            if length(record.sequence) != length(record.quality)
                throw(ArgumentError("the length of sequence does not match the length of quality"))
            end
            BioCore.ReaderHelper.resize_and_copy!(record.data, data, p-(@pos)+1:p-1)
            if !isempty(header2_identifier)
                check_identical(record.data, record.identifier, record.data, header2_identifier)
            end
            if !isempty(header2_description)
                check_identical(record.data, record.description, record.data, header2_description)
            end
            record.filled = 1:(@pos)-1
            state.read = true
            @unmark
            if p ≤ p_end && data[p] == UInt8('@')
                # cancel lookahead
                p -= 1
            end
            @escape
        end)

    # Generate readrecord! function.
    eval(
        BioCore.ReaderHelper.generate_readrecord_function(
            Record, machine, actions, initcode, exitcode, loopunroll=10))
end
