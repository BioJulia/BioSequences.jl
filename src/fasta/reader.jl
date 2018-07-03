# FASTA Reader
# ============

mutable struct Reader <: BioCore.IO.AbstractReader
    state::BioCore.Ragel.State
    index::Union{Index, Nothing}

    function Reader(input::BufferedStreams.BufferedInputStream, index)
        return new(BioCore.Ragel.State(file_machine.start_state, input), index)
    end
end

"""
    FASTA.Reader(input::IO; index = nothing)

Create a data reader of the FASTA file format.

# Arguments
* `input`: data source
* `index=nothing`: filepath to a random access index (currently *fai* is supported)
"""
function Reader(input::IO; index = nothing)
    if isa(index, AbstractString)
        index = Index(index)
    else
        if index != nothing
            throw(ArgumentError("index must be a filepath or nothing"))
        end
    end
    return Reader(BufferedStreams.BufferedInputStream(input), index)
end

function Base.eltype(::Type{Reader})
    return Record
end

function BioCore.IO.stream(reader::Reader)
    return reader.state.stream
end

function Base.getindex(reader::Reader, name::AbstractString)
    if reader.index == nothing
        throw(ArgumentError("no index attached"))
    end
    seekrecord(reader.state.stream, reader.index, name)
    reader.state.cs = file_machine.start_state
    reader.state.finished = false
    return read(reader)
end

isinteractive() && @info "Compiling FASTA parser..."
const record_machine, file_machine = (function ()
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
    identifier.actions[:enter] = [:mark]
    identifier.actions[:exit]  = [:identifier]

    description = cat(any() \ hspace, re"[^\r\n]*")
    description.actions[:enter] = [:mark]
    description.actions[:exit]  = [:description]

    header = cat('>', identifier, opt(cat(rep1(hspace), description)))
    header.actions[:enter] = [:anchor]
    header.actions[:exit]  = [:header]

    # '*': terminal, `-': gap
    letters = re"[A-Za-z*\-]+"
    letters.actions[:enter] = [:movable_anchor]
    letters.actions[:exit]  = [:letters]

    sequence = opt(cat(letters, rep(cat(rep1(alt(hspace, newline)), letters))))
    sequence.actions[:enter] = [:sequence_start]

    record = cat(header, rep1(newline), sequence, rep1(newline))
    record.actions[:exit] = [:record]

    record_trailing = cat(header, rep1(newline), sequence)
    record_trailing.actions[:exit] = [:record]

    file = cat(rep(newline), rep(record), opt(record_trailing))

    return map(Automa.compile, (record, file))
end)()

#=
write("fasta.dot", Automa.machine2dot(file_machine))
run(`dot -Tsvg -o fasta.svg fasta.dot`)
=#

const record_actions = Dict(
    :identifier  => :(record.identifier  = (mark:p-1)),
    :description => :(record.description = (mark:p-1)),
    :header => quote
        @assert record.data === data
        copyto!(record.data, 1, record.data, 1, p - 1)
        filled += p - 1
        if p ≤ lastindex(data)
            data[p] = UInt8('\n')
            filled += 1
        end
    end,
    :letters => quote
        let len = p - mark
            copyto!(record.data, filled + 1, record.data, mark, p - mark)
            filled += len
            record.sequence = first(record.sequence):last(record.sequence)+len
        end
    end,
    :sequence_start => :(record.sequence = filled+1:filled),
    :record => :(record.filled = 1:filled),
    :anchor => :(),
    :movable_anchor => :(mark = p),
    :mark => :(mark = p),
    :countline => :(#= linenum += 1 =#))
eval(
    BioCore.ReaderHelper.generate_index_function(
        Record,
        record_machine,
        :(filled = mark = 0),
        record_actions))
eval(
    BioCore.ReaderHelper.generate_read_function(
        Reader,
        file_machine,
        :(filled = mark = 0),
        merge(record_actions, Dict(
            :identifier  => :(record.identifier  = (mark:p-1) .- stream.anchor .+ 1),
            :description => :(record.description = (mark:p-1) .- stream.anchor .+ 1),
            :header => quote
                range = BioCore.ReaderHelper.upanchor!(stream):p-1
                BioCore.ReaderHelper.resize_and_copy!(record.data, filled + 1, data, range)
                filled += length(range)
                if filled + 1 ≤ lastindex(record.data)
                    record.data[filled+1] = UInt8('\n')
                else
                    push!(record.data, UInt8('\n'))
                end
                filled += 1
            end,
            :sequence_start => :(record.sequence = filled+1:filled),
            :letters => quote
                let len = p - stream.anchor
                    BioCore.ReaderHelper.append_from_anchor!(record.data, filled + 1, stream, p - 1)
                    filled += len
                    record.sequence = first(record.sequence):last(record.sequence)+len
                end
            end,
            :record => quote
                record.filled = 1:filled
                found_record = true
                @escape
            end,
            :countline => :(linenum += 1),
            :anchor => :(BioCore.ReaderHelper.anchor!(stream, p)),
            :movable_anchor => :(BioCore.ReaderHelper.anchor!(stream, p, false))))))
