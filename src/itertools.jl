# TODO: This code should be placed somewhere else (BioCore.jl?, maybe).

module ReaderIterTools

abstract type AbstractRecordIterator{T} end

function eachrecord end

mutable struct RecordIteratorState{T}
    # machine state
    state::Int
    # line number
    linenum::Int
    # is record filled?
    filled::Bool
    # placeholder
    record::T
end

function Base.iteratorsize(::Type{<:AbstractRecordIterator})
    return Base.SizeUnknown()
end

function Base.eltype(::Type{AbstractRecordIterator{T}}) where T
    return T
end

function Base.start(iter::AbstractRecordIterator{T}) where T
    return ReaderIterTools.RecordIteratorState(1, 1, false, T())
end

function Base.next(iter::AbstractRecordIterator{T}, state) where T
    @assert state.filled
    record = copy(state.record)
    state.filled = false
    return record, state
end


end
