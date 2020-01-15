# Specialized printing/showing methods
#

function Base.print(io::IO, seq::LongSequence{A}; width::Integer = 0) where {A<:Alphabet}
    return _print(io, seq, width, codetype(A()))
end

# Dispatch to generic method in biosequences/printing.jl
function _print(io::IO, seq::LongSequence{<:Alphabet}, width::Integer, ::AlphabetCode)
    return _print(io, seq, width)
end

# Specialized method for ASCII alphabet
function _print(io::IO, seq::LongSequence{<:Alphabet}, width::Integer, ::AsciiAlphabet)
    # I don't like to have to do this, but in Julia 1.3, system buffers are IO-locked.
    if (width < 1) | (length(seq) ≤ width)
        return print(io, String(seq))
    end
    if length(seq) < 4096
        buffer = SimpleBuffer(io, padded_length(length(seq), width))
    else
        buffer = SimpleBuffer(io)
    end
    return _print(buffer, seq, width, AsciiAlphabet())
end

function _print(buffer::SimpleBuffer, seq::LongSequence{<:Alphabet}, width::Integer, ::AsciiAlphabet)
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
