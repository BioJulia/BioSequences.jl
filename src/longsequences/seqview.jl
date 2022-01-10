################################### Construction etc

# Mention in docs no boundscheck is performed on instantiation
"""
	LongSubSeq{A <: Alphabet}

A view into a `LongSequence`. This shares data buffer with the underlying
sequence, and is therefore much faster to instantiate than a `LongSequence`.
Modifying the view changes the sequence and vice versa.

# Examples
```jldoctest
julia> LongSubSeq(dna"TAGA", 2:3)
AG

julia view(dna"TAGA", 2:3)
AG
```
"""
struct LongSubSeq{A<:Alphabet} <: BioSequence{A}
    data::Vector{UInt64}
    part::UnitRange{Int}

	# Added to reduce method ambiguities
	LongSubSeq{A}(data::Vector{UInt64}, part::UnitRange{Int}) where A = new{A}(data, part)
end

# These unions are significant because LongSubSeq and LongSequence have the same
# encoding underneath, so many methods can be shared.
const SeqOrView{A} = Union{LongSequence{A}, LongSubSeq{A}}
const NucleicSeqOrView = SeqOrView{<:NucleicAcidAlphabet}

Base.length(v::LongSubSeq) = last(v.part) - first(v.part) + 1
Base.copy(v::LongSubSeq{A}) where A = LongSequence{A}(v)

encoded_data_eltype(::Type{<:LongSubSeq}) = encoded_data_eltype(LongSequence)
symbols_per_data_element(x::LongSubSeq) = div(64, bits_per_symbol(Alphabet(x)))

@inline function bitindex(x::LongSubSeq, i::Integer)
    N = BitsPerSymbol(Alphabet(typeof(x)))
    bitindex(N, encoded_data_eltype(typeof(x)), i % UInt + first(x.part) - 1)
end

# Constructors
function LongSubSeq{A}(seq::LongSequence{A}) where A
    return LongSubSeq{A}(seq.data, 1:length(seq))
end

function LongSubSeq{A}(seq::LongSubSeq{A}) where A
    return LongSubSeq{A}(seq.data, seq.part)
end

function LongSubSeq{A}(seq::LongSequence{A}, part::AbstractUnitRange{<:Integer}) where A
    @boundscheck checkbounds(seq, part)
    return LongSubSeq{A}(seq.data, UnitRange{Int}(part))
end

function LongSubSeq{A}(seq::LongSubSeq{A}, part::AbstractUnitRange{<:Integer}) where A
    @boundscheck checkbounds(seq, part)
    newpart = first(part) + first(seq.part) - 1 : last(part) + first(seq.part) - 1
    return LongSubSeq{A}(seq.data, newpart)
end

function LongSubSeq(seq::SeqOrView{A}, i) where A
	return LongSubSeq{A}(seq, i)
end

LongSubSeq(seq::SeqOrView, ::Colon) = LongSubSeq(seq, 1:lastindex(seq))
LongSubSeq(seq::BioSequence{A}) where A = LongSubSeq{A}(seq)

Base.view(seq::SeqOrView, part::AbstractUnitRange) = LongSubSeq(seq, part)

# Conversion
function LongSequence(s::LongSubSeq{A}) where A
	_copy_seqview(LongSequence{A}, s)
end

function LongSequence{A}(seq::LongSubSeq{A}) where {A<:NucleicAcidAlphabet}
	_copy_seqview(LongSequence{A}, seq)
end

function _copy_seqview(T, s::LongSubSeq)
	first = firstbitindex(s)
	v = s.data[index(first):index(lastbitindex(s))]
	res = T(v, length(s) % UInt)
	return zero_offset!(res, offset(first) % UInt)
end

function (::Type{T})(seq::LongSequence{<:NucleicAcidAlphabet{N}}) where
	{N, T<:LongSubSeq{<:NucleicAcidAlphabet{N}}}
	T(seq.data, 1:length(seq))
end

function Base.convert(::Type{T1}, seq::T2) where
	{T1 <: Union{LongSequence, LongSubSeq}, T2 <: Union{LongSequence, LongSubSeq}}
	return T1(seq)
end

# Indexing
function Base.getindex(seq::LongSubSeq, part::AbstractUnitRange{<:Integer})
	return LongSubSeq(seq, part)
end

Base.parentindices(seq::LongSubSeq) = (seq.part,)
