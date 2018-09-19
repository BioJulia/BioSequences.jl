# Printing, show and parse
# ------------------------

function Base.print(io::IO, seq::BioSequence; width::Integer = 0)
    col = 1
    for x in seq
        if width > 0 && col > width
            write(io, '\n')
            col = 1
        end
        print(io, x)
        col += 1
    end
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
            for x in seq
                print(io, convert(Char, x))
            end
        end
    end
end

function string_compact(seq::BioSequence)
    buf = IOBuffer()
    showcompact(buf, seq)
    return String(take!(buf))
end

Base.parse(::Type{S}, str::AbstractString) where {S<:BioSequence} = convert(S, str)