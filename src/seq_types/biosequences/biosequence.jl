# BioSequence
# ===========
#
# Abstract biological sequence type.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md


# Type Definition
# ===============

"""
# An abstract biological sequence type.

## Required methods, traits & interfaces

Any subtype `S <: BioSequence` should implement the following methods:

* `alphabet_t(seq::S)`: return the type of the alphabet of `seq`.
* `Base.length(seq::S)`: return the length of `seq`.
* `inbounds_getindex(seq::S, i::Integer)`: return the element at `i` of `seq`
  without checking bounds.
* `find_next_ambiguous(seq::S)`: return the index of the next ambiguous symbol in the
  sequence.

## Provided, methods traits & interfaces

As a result you can expect any variable of type S <: BioSequence to implement
the above `Base` methods, as well as the following `Base` methods which are
implemetned in BioSequences.jl for any S<:BioSequence which implements said
required `Base` methods above.

* `Base.eltype(::Type{S})`: return the element type of `S`.
* `Base.size(seq::S)`: Return the size of `seq`.
* `Base.endof(seq::S)`: Return the last index of `seq`.
* `Base.eachindex(seq::S)`: Return an iterator over all indicies of `seq`.
* `Base.getindex(seq::S, i::Integer)`: Return the biological symbol `seq`
  contains at position `i`.
* `Base.:(==)(seq1::S, seq2::S)`: Check if `seq1` and `seq2` can be considered
  equal in value.
* `Base.isless(seq1::S, seq2::S)`: Check if `seq1` is considered lesser in
  value than `seq2`.
* `Base.cmp(seq1::S, seq2::S)`: Returns -1, 0, or 1 depending on whether
  `seq1` is less than, equal to, or greater than `seq2`, respectively.
* `Base.findnext(seq::S, val, start::Integer)`: Return the index of the next
  occurance of `val` in `seq`, starting from index `start`.
* `Base.findprev(seq::S, val, start::Integer)`: Return the index of the
  previous occurance of `val` in `seq`, starting from index `start`.
* `Base.findfirst(seq::S, val)`: Return the index of the first occurance of
  `val` in `seq`.
* `Base.findlast(seq::S, val)`: Return the index of the final occurance of
  `val` in `seq`.
* `Base.isempty(seq::S)`: Determine whether `seq` has no elements.
"""

abstract type BioSequence end

# This is useful for obscure reasons. We use SeqRecord{BioSequence} for reading
# sequence in an undetermined alphabet, but a consequence that we need to be
# able to construct a `Sequence`.
function BioSequence()
    return DNASequence()
end


# Base Methods
# ============


# Indexing and iteration
# ----------------------

Base.eltype(::Type{T}) where T <: BioSequence = eltype(alphabet_t(T))
Base.eltype(seq::BioSequence) = eltype(alphabet_t(seq))
Base.size(seq::BioSequence) = (length(seq),)
Base.lastindex(seq::BioSequence) = length(seq)
Base.eachindex(seq::BioSequence) = 1:lastindex(seq)

@inline function Base.checkbounds(seq::BioSequence, i::Integer)
    if 1 ≤ i ≤ endof(seq)
        return true
    end
    throw(BoundsError(seq, i))
end

"""
Return the `Alpahbet` type defining the possible biological symbols
and their encoding for a given biological sequence.
"""
@inline function alphabet_t(::Type{S}) where S <: BioSequence
    error(string("This sequence type trait has not been defined for BioSequence type: ", S))
end

# This version of alphabet_t is automatically defined for any BioSequence type, is more for conveinience.
@inline function alphabet_t(seq::BioSequence)
    return alphabet_t(typeof(seq))
end

@inline function bitindex_t(seq::BioSequence)
    error(string("This sequence type trait has not been defined for BioSequence type: ", typeof(seq)))
end

"""
Return the data member of `seq` that stores the encoded sequence data.
"""
@inline function encoded_data(seq::BioSequence)
    error(
        string(
            "encoded_data has not been defined for BioSequence type: ",
            typeof(seq),
            ". It is required for any BioSequence subtype."
        )
    )
end

# Bit indexing
# ------------

# Bit indexing biosequence traits and trait-like methods...
bits_per_symbol(::Type{S}) where {S <: BioSequence} = bits_per_symbol(alphabet_t(S))
bits_per_symbol(seq::BioSequence) = bits_per_symbol(alphabet_t(seq))

bits_per_symbol_t(seq::BioSequence) = bits_per_symbol_t(alphabet_t(seq))

encoded_data_eltype(seq::BioSequence) = eltype(encoded_data(seq))

@inline function symbols_per_data_element(seq::BioSequence)
    return div(8 * sizeof(encoded_data_eltype(seq)), bits_per_symbol(seq))
end

@inline function bitindex(seq::BioSequence, i::Integer)
    return bitindex(bits_per_symbol_t(seq), encoded_data_eltype(seq), i)
end

@inline function bindata_mask(seq::BioSequence)
    return bitmask(alphabet_t(seq))
end








# Base Methods
# ============

# Indexing and iteration
# ----------------------

include("indexing.jl")

function Base.iterate(seq::BioSequence, i::Int = 1)
    if i > lastindex(seq)
        return nothing
    else
        return inbounds_getindex(seq, i), i + 1
    end
end


# Predicates & comparisons
# ------------------------

function Base.cmp(seq1::BioSequence, seq2::BioSequence)
    m = endof(seq1)
    n = endof(seq2)
    for i in 1:min(m, n)
        c = cmp(inbounds_getindex(seq1, i),
                inbounds_getindex(seq2, i))
        if c != 0
            return c
        end
    end
    return cmp(m, n)
end

function Base.:(==)(seq1::BioSequence, seq2::BioSequence)
    return eltype(seq1)    == eltype(seq2) &&
           length(seq1)    == length(seq2) &&
           cmp(seq1, seq2) == 0
end

Base.isless(seq1::BioSequence, seq2::BioSequence) = cmp(seq1, seq2) < 0

function Base.isempty(seq::BioSequence)
    return length(seq) == 0
end


# Finders
# -------

function Base.findnext(val, seq::BioSequence, start::Integer)
    checkbounds(seq, start)
    v = convert(eltype(seq), val)
    for i in start:lastindex(seq)
        if inbounds_getindex(seq, i) == v
            return i
        end
    end
    return nothing
end

function Base.findprev(val, seq::BioSequence, start::Integer)
    checkbounds(seq, start)
    v = convert(eltype(seq), val)
    for i in start:-1:1
        if inbounds_getindex(seq, i) == v
            return i
        end
    end
    return nothing
end

Base.findfirst(val, seq::BioSequence) = findnext(val, seq, 1)
Base.findlast(val, seq::BioSequence) = findprev(val, seq, lastindex(seq))

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

include("operations.jl")
