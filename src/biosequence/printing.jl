# Printing, show and parse
# ------------------------

Base.summary(seq::BioSequence{<:DNAAlphabet}) = string(length(seq), "nt ", "DNA Sequence")
Base.summary(seq::BioSequence{<:RNAAlphabet}) = string(length(seq), "nt ", "RNA Sequence")
Base.summary(seq::BioSequence{<:AminoAcidAlphabet}) = string(length(seq), "aa ", "Amino Acid Sequence")

# Buffer type. Not exposed to user, so code should be kept simple and performant.
# B is true if it is buffered and false if it is not
mutable struct SimpleBuffer{B, T} <: IO
    len::UInt
    arr::Vector{UInt8}
    io::T
end

SimpleBuffer(io::IO) = SimpleBuffer{true, typeof(io)}(0, Vector{UInt8}(undef, 1024), io)
SimpleBuffer(io::IO, len) = SimpleBuffer{false, typeof(io)}(0, Vector{UInt8}(undef, len), io)

function Base.write(sb::SimpleBuffer{true}, byte::UInt8)
    sb.len ≥ 1024 && flush(sb)
    sb.len += 1
    @inbounds sb.arr[sb.len] = byte
end

function Base.write(sb::SimpleBuffer{false}, byte::UInt8)
    len = sb.len + 1
    sb.len = len
    @inbounds sb.arr[len] = byte
end

# Flush entire buffer to its io
@noinline function Base.flush(sb::SimpleBuffer{true})
    arr = sb.arr
    GC.@preserve arr unsafe_write(sb.io, pointer(arr), UInt(1024))
    sb.len = 0
end

# Flush all unflushed bytes to its io. Does not close source or make buffer unwritable
function Base.close(sb::SimpleBuffer)
    iszero(sb.len) && return 0
    arr = sb.arr
    GC.@preserve arr unsafe_write(sb.io, pointer(arr), sb.len)
    sb.len = 0
end

function padded_length(len::Integer, width::Integer)
    den = ifelse(width < 1, typemax(Int), width)
    return len + div(len-1, den)
end

function Base.print(io::IO, seq::BioSequence{A}; width::Integer = 0) where {A<:Alphabet}
    return _print(io, seq, width, codetype(A()))
end

# Generic method. The different name allows subtypes of BioSequence to
# selectively call the generic print despite being more specific type
function _print(io::IO, seq::BioSequence, width::Integer, ::UnicodeAlphabet)
    col = 0
    for x in seq
        col += 1
        if (width > 0) & (col > width)
            write(io, '\n')
            col = 1
        end
        print(io, x)
    end
    return nothing
end

# Specialized method for ASCII alphabet
function _print(io::IO, seq::BioSequence, width::Integer, ::AsciiAlphabet)
    # If seq is large, always buffer for memory efficiency
    if length(seq) ≥ 4096
        return _print(SimpleBuffer(io), seq, width, AsciiAlphabet())
    end
    if (width < 1) | (length(seq) ≤ width)
        # Fastest option
        return print(io, String(seq))
    else
        # Second fastest option
        buffer = SimpleBuffer(io, padded_length(length(seq), width))
        return _print(buffer, seq, width, AsciiAlphabet())
    end
end

function _print(buffer::SimpleBuffer, seq::BioSequence, width::Integer, ::AsciiAlphabet)
    col = 0
    @inbounds for i in eachindex(seq)
        col += 1
        if (width > 0) & (col > width)
            write(buffer, UInt8('\n'))
            col = 1
        end
        write(buffer, stringbyte(seq[i]))
    end
    close(buffer)
    return nothing
end

Base.show(io::IO, seq::BioSequence) = showcompact(io, seq)

function Base.show(io::IO, ::MIME"text/plain", seq::BioSequence)
    println(io, summary(seq), ':')
    showcompact(io, seq)
end

function showcompact(io::IO, seq::BioSequence)
    # don't show more than this many characters
    # to avoid filling the screen with junk
    if isempty(seq)
        print(io, "< EMPTY SEQUENCE >")
    else
        width = displaysize()[2]
        if length(seq) > width
            half = div(width, 2)
            for i in 1:half-1
                print(io, seq[i])
            end
            print(io, '…')
            for i in lastindex(seq)-half+2:lastindex(seq)
                print(io, seq[i])
            end
        else
            print(io, seq)
        end
    end
end

function string_compact(seq::BioSequence)
    buf = IOBuffer()
    showcompact(buf, seq)
    return String(take!(buf))
end
