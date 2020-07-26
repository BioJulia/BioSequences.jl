################################### Construction etc

# Mention in docs no boundscheck is performed on instantiation
struct SeqView{A<:Alphabet} <: BioSequence{A}
    data::Vector{UInt64}
    part::UnitRange{Int}
end

# These unions are significant because SeqView and LongSequence have the same
# encoding underneath, so many methods can be shared.
const SeqOrView{A} = Union{LongSequence{A}, SeqView{A}}
const NucleicSeqOrView = SeqOrView{<:NucleicAcidAlphabet}

Base.length(v::SeqView) = last(v.part) - first(v.part) + 1
encoded_data(v::SeqView) = v.data

function SeqView(seq::LongSequence{A}, part::UnitRange{Int}) where A
    @boundscheck checkbounds(seq, part)
    return SeqView{A}(seq.data, part)
end

Base.view(seq::LongSequence, part::UnitRange) = SeqView(seq, part)
function Base.view(seq::SeqView, part::UnitRange)
	offset = first(seq.part) - 1
	return SeqView(seq, first(part)+offset : last(part) + offset)
end

function LongSequence(s::SeqView{A}) where A
	_copy_seqview(LongSequence{A}, s)
end

function (::Type{T})(seq::SeqView{<:NucleicAcidAlphabet{N}}) where
	{N, T<:LongSequence{<:NucleicAcidAlphabet{N}}}
	_copy_seqview(T, seq)
end

function _copy_seqview(T, seq::SeqView)
	first = firstbitindex(s)
	v = s.data[index(first):index(lastbitindex(s))]
	rightshift!(v, offset(first))
	return T(v, length(s))
end

function (::Type{T})(seq::LongSequence{<:NucleicAcidAlphabet{N}}) where
	{N, T<:SeqView{<:NucleicAcidAlphabet{N}}}
	T(seq.data, 1:length(seq))
end

function Base.convert(::Type{T1}, seq::T2) where
	{T1 <: Union{LongSequence, SeqView}, T2 <: Union{LongSequence, SeqView}}
	return T1(seq)
end
