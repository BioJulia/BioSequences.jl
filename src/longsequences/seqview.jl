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

# Constructors
function SeqView{A}(seq::LongSequence{A}) where {A <: Alphabet}
    return SeqView{A}(seq.data, 1:length(seq))
end

function SeqView{A}(seq::SeqView{A}) where {A <: Alphabet}
    return SeqView{A}(seq.data, seq.part)
end

function SeqView(seq::SeqOrView{A}) where {A <: Alphabet}
	return SeqView{A}(seq)
end

function SeqView{A}(seq::LongSequence{A}, part::UnitRange{Int}) where {A <: Alphabet}
    @boundscheck checkbounds(seq, part)
    return SeqView{A}(seq.data, part)
end

function SeqView{A}(seq::SeqView{A}, part::UnitRange{<:Integer}) where {A <: Alphabet}
    @boundscheck checkbounds(seq, part)
    newpart = first(part) + first(seq.part) - 1 : last(part) + first(seq.part) - 1
    return SeqView{A}(seq.data, newpart)
end

function SeqView(seq::SeqOrView{A}, i) where {A <: Alphabet}
	return SeqView{A}(seq, i)
end

Base.view(seq::SeqOrView, part::UnitRange) = SeqView(seq, part)

# Conversion
function LongSequence(s::SeqView{A}) where {A <: Alphabet}
	_copy_seqview(LongSequence{A}, s)
end

function (::Type{T})(seq::SeqView{<:NucleicAcidAlphabet{N}}) where
	{N, T<:LongSequence{<:NucleicAcidAlphabet{N}}}
	_copy_seqview(T, seq)
end

function _copy_seqview(T, s::SeqView)
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


# Indecing
function Base.getindex(seq::SeqView, part::UnitRange{<:Integer})
	return SeqView(seq, part)
end
