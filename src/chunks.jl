struct ChunkIterator{S}
    v::Vector{UInt64}
    last_chunk_index::Int
end

function Base.iterate(it::ChunkIterator{<:LongSequence}, state::Int=1)
    state > it.last_chunk_index && return nothing
    (@inbounds(it.v[state]), state + 1)
end

Base.length(it::ChunkIterator{<:LongSequence}) = it.last_chunk_index

function iter_chunks(seq::S)::Tuple{ChunkIterator{S}, Tuple{UInt64, UInt8}} where {S <: LongSequence}
    v = seq.data
    @assert length(v) == cld(length(seq), symbols_per_data_element(seq))
    remainder = (UInt(length(seq)) % UInt(64)) % UInt8
    start = if iszero(remainder)
        (UInt64(0), 0x00)
    else
        (@inbounds(v[end]), remainder)
    end
    (
        ChunkIterator{typeof(seq)}(v, lastindex(v) - !iszero(remainder)),
        start
    )
end
