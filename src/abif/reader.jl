# ABIF Reader
# ===========


export get_tags, tagelements

struct AbifDirEntry
    name::String
    number::Int32
    element_type::Int32
    element_size::Int32
    num_elements::Int32
    data_size::Int32
    data_offset::Int32
end

"""
    ABIF.Reader(input::IO)
Create a data reader of the ABIF file format.
# Arguments
* `input`: data source
"""
mutable struct Reader{T<:IO} <: BioCore.IO.AbstractReader
    # input stream
    input::T

    # Tags
    dirs::Vector{AbifDirEntry}
end

function Reader(input::IO)
    header = read_abif_header(input)
    tags   = read_abif_tags(input, header)

    return Reader(input, tags)
end

function BioCore.IO.stream(reader::Reader)
    return reader.input
end

function Base.iterate(p::Reader, k::Int=1)
    if k > length(p)
        return nothing
    else
        return getindex(p, p.dirs[k]), k + 1
    end
end

function Base.length(a::Reader)
    return length(a.dirs)
end

function Base.checkbounds(p::Reader, k::Integer)
    if 1 ≤ k ≤ length(p)
        return true
    end
    throw(BoundsError(p, k))
end

function Base.getindex(a::Reader, k::Integer)
    checkbounds(a, k)
    return Dict([parse_data_tag(a, a.dirs[k])])
end

function Base.getindex(a::Reader, t::AbstractString)
    tag = filter((x) -> x.name == t, a.dirs)

    if isempty(tag)
        throw(KeyError(t))
    end

    if length(tag) > 1
        return Dict([parse_data_tag(a, b) for b in tag])
    end

    return Dict([parse_data_tag(a, first(tag))])
end

function Base.getindex(a::Reader, t::AbifDirEntry)
    return Dict([parse_data_tag(a, t)])
end

function Base.getindex(a::Reader, t::Array{AbifDirEntry})
    return Dict([parse_data_tag(a, k) for k in t])
end

"""
    get_tags(input::Reader)
Returns all existing tags
# Arguments
* `input`: Reader
"""
function get_tags(a::Reader)
    return [tag for tag in a.dirs]
end

"""
    get_tags(input::Reader, tag_name::AbstractString)
Returns all existing tags by name
# Arguments
* `input`: Reader
* `tag_name`: AbstractString
"""
function get_tags(a::Reader, t::AbstractString)
    return [tag for tag in a.dirs if isequal(tag.name, t)]
end

"""
    tagelements(stream::Reader, tag_name::AbstractString)
Returns the number of how many tags exists by the same name
# Arguments
* `stream`: Reader
* `tag_name`: AbstractString
"""
function tagelements(a::Reader, t::AbstractString)
    return length(get_tags(a, t))
end

# extract the header
function read_abif_header(input::IO)
    seekstart(input)
    signature = read(input, 4)

    if is_abif_signature(signature)
        version = ntoh(first(read(input, Int16)))
        return parse_directory(input, position(input))
    else
        error("Invalid File Signature")
    end
end

# extract all tags in file
function read_abif_tags(input::IO, header::AbifDirEntry)
    tags = AbifDirEntry[]
    for index in 0:header.num_elements-1
        start = header.data_offset + index * header.element_size
        tag = parse_directory(input, start)
        push!(tags, tag)
    end
    return tags
end

# extract DirEtry
function parse_directory(input::IO, pos::Int64)
    seek(input, pos)

    name         = String(read(input, 4))
    number       = ntoh(read(input, Int32))
    element_type = ntoh(read(input, UInt16))
    element_size = ntoh(read(input, Int16))
    num_elements = ntoh(read(input, Int32))
    data_size    = ntoh(read(input, Int32))

    if data_size <= 4
        data_offset = position(input)
    else
        data_offset = ntoh(read(input, Int32))
    end

    AbifDirEntry(name, number, element_type, element_size, num_elements, data_size, data_offset)
end

# read bytes according to the element type, other values are unsupported or legacy.
function parse_data_tag(a::Reader, tag::AbifDirEntry)
    if tag.data_offset > 0
        seek(a.input, tag.data_offset)

        data = Tuple{String, Union{AbstractString, Integer}}

        if tag.element_type == 2
            data = String(read(a.input, tag.data_size))

        elseif tag.element_type == 4
            data = [read(a.input, Int16) for _ in 1:tag.num_elements]
            data = convert_to_int(data)

        elseif tag.element_type == 5
            data = [read(a.input, Int32) for _ in 1:tag.num_elements]
            data = convert_to_int(data)

        elseif tag.element_type == 7
            data = [read(a.input, Float32) for _ in 1:tag.num_elements]
            data = convert_to_float(data)

        elseif tag.element_type == 10
            year   = ntoh(read(a.input, Int16))
            month  = Int(read(a.input, UInt8))
            day    = Int(read(a.input, UInt8))
            data   = "$year-$month-$day"

        elseif tag.element_type == 11
            hour    = Int(ntoh(read(a.input, UInt8)))
            minute  = Int(ntoh(read(a.input, UInt8)))
            second  = Int(ntoh(read(a.input, UInt8)))
            hsecond = Int(ntoh(read(a.input, UInt8)))
            data    = "$hour:$minute:$second:$hsecond"

        elseif tag.element_type == 13
            data = Bool(first(read(a.input, tag.data_size)))

        elseif tag.element_type == 18
            data = String(read(a.input, tag.data_size)[2:end])

        elseif tag.element_type == 19
            data = String(read(a.input, tag.data_size)[1:end-1])

        else
            data = read(a.input, tag.data_size)
        end

        return format_data!(data, tag, tagelements(a, tag.name))
    end
end

# if TAG has more than one element, concatenate element number on name
function format_data!(data::Any, tag::AbifDirEntry, elements::Integer)
    if elements > 1
        return ("$(tag.name)$(tag.number)", data)
    end
    return ("$(tag.name)", data)
end

# convert a array of big-endian values
function convert_to_int(data::AbstractArray)
    result = Int[]

    for d in data
        r = Int(ntoh(d))
        push!(result, r)
    end

    if length(result) <= 1
        return first(result)
    end

    return result
end

# convert a array of big-endian values
function convert_to_float(data::AbstractArray)
    result = Float32[]
    for d in data
        r = Float32(ntoh(d))
        push!(result, r)
    end

     if length(result) <= 1
        return first(result)
     end

    return result
end

# Check if the file has a valid Signature
function is_abif_signature(signature::Array{UInt8})
    if signature != b"ABIF"
        return false
    end
    return true
end
