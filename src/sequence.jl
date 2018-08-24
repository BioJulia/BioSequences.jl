# Sequence
# ========
#
# Abstract biological sequence type.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

"""
Abstract biological sequence type.

Any subtype `S <: Sequence` should implement the following methods:

* `Base.length(seq::S)`: return the length of `seq`.
* `Base.eltype(::Type{S})`: return the element type of `S`.
* `inbounds_getindex(seq::S, i::Integer)`: return the element at `i` of `seq`
    without checking bounds
"""
abstract type Sequence end

# This is useful for obscure reasons. We use SeqRecord{Sequence} for reading
# sequence in an undetermined alphabet, but a consequence that we need to be
# able to construct a `Sequence`.
function Sequence()
    return DNASequence()
end

Base.size(seq::Sequence) = (length(seq),)
Base.lastindex(seq::Sequence) = length(seq)
Base.eachindex(seq::Sequence) = 1:lastindex(seq)

@inline function Base.checkbounds(seq::Sequence, i::Integer)
    if 1 ≤ i ≤ lastindex(seq)
        return true
    end
    throw(BoundsError(seq, i))
end

@inline function Base.getindex(seq::Sequence, i::Integer)
    @boundscheck checkbounds(seq, i)
    return inbounds_getindex(seq, i)
end

function Base.iterate(seq::Sequence, i::Int=1)
    if i > lastindex(seq)
        return nothing
    else
        return inbounds_getindex(seq, i), i + 1
    end
end


# Comparison
# ----------

function Base.:(==)(seq1::Sequence, seq2::Sequence)
    return eltype(seq1)    == eltype(seq2) &&
           length(seq1)    == length(seq2) &&
           cmp(seq1, seq2) == 0
end

Base.isless(seq1::Sequence, seq2::Sequence) = cmp(seq1, seq2) < 0

function Base.cmp(seq1::Sequence, seq2::Sequence)
    m = length(seq1)
    n = length(seq2)
    for i in 1:min(m, n)
        c = cmp(inbounds_getindex(seq1, i),
                inbounds_getindex(seq2, i))
        if c != 0
            return c
        end
    end
    return cmp(m, n)
end


# Finders
# -------

function Base.findnext(val, seq::Sequence, start::Integer)
    checkbounds(seq, start)
    v = convert(eltype(seq), val)
    for i in start:lastindex(seq)
        if inbounds_getindex(seq, i) == v
            return i
        end
    end
    return nothing
end

function Base.findprev(val, seq::Sequence, start::Integer)
    checkbounds(seq, start)
    v = convert(eltype(seq), val)
    for i in start:-1:1
        if inbounds_getindex(seq, i) == v
            return i
        end
    end
    return nothing
end

Base.findfirst(val, seq::Sequence) = findnext(val, seq, 1)
Base.findlast(val, seq::Sequence) = findprev(val, seq, lastindex(seq))


# GC content
# ----------

"""
    gc_content(seq::Sequence)

Calculate GC content of `seq`.
"""
function gc_content(seq::Sequence)
    if !(eltype(seq) <: NucleicAcid)
        throw(ArgumentError("not a nucleic acid sequence"))
    end
    if isempty(seq)
        return 0.0
    else
        return count_gc(seq) / length(seq)
    end
end

function count_gc(seq::Sequence)
    return count(isGC, seq)
end


# Predicates
# ----------

function Base.isempty(seq::Sequence)
    return length(seq) == 0
end

"""
    ispalindromic(seq::Sequence)

Return `true` if `seq` is a palindromic sequence; otherwise return `false`.
"""
function ispalindromic(seq::Sequence)
    if !(eltype(seq) <: NucleicAcid)
        error("elements must be nucleotide")
    end

    for i in 1:cld(length(seq), 2)
        if seq[i] != complement(seq[end-i+1])
            return false
        end
    end

    return true
end

"""
    hasambiguity(seq::Sequence)

Return `true` if `seq` has an ambiguous symbol; otherwise return `false`.
"""
function hasambiguity(seq::Sequence)
    for x in seq
        if isambiguous(x)
            return true
        end
    end
    return false
end

"""
    isrepetitive(seq::Sequence, n::Integer=length(seq))

Return `true` if and only if `seq` contains a repetitive subsequence of length `≥ n`.
"""
function isrepetitive(seq::Sequence, n::Integer=length(seq))
    if n < 0
        error("repetition must be non-negative")
    elseif isempty(seq)
        return n == 0
    end

    rep = 1
    if rep ≥ n
        return true
    end
    last = first(seq)
    for i in 2:lastindex(seq)
        x = seq[i]
        if x == last
            rep += 1
            if rep ≥ n
                return true
            end
        else
            rep = 1
        end
        last = x
    end

    return false
end


# Transformations
# ---------------

"Create a copy of a sequence with gap characters removed."
ungap(seq::Sequence)  =  filter(x -> x != BioSymbols.gap(eltype(seq)), seq)

"Remove gap characters from a sequence. Modifies the input sequence."
ungap!(seq::Sequence) = filter!(x -> x != BioSymbols.gap(eltype(seq)), seq)


# Printers
# --------

function Base.print(io::IO, seq::Sequence; width::Integer = 0)
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

Base.show(io::IO, seq::Sequence) = showcompact(io, seq)

function Base.show(io::IO, ::MIME"text/plain", seq::Sequence)
    println(io, summary(seq), ':')
    showcompact(io, seq)
end

function showcompact(io::IO, seq::Sequence)
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

function string_compact(seq::Sequence)
    buf = IOBuffer()
    showcompact(buf, seq)
    return String(take!(buf))
end

Base.parse(::Type{S}, str::AbstractString) where {S<:Sequence} = convert(S, str)
