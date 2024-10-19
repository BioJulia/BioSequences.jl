struct LongSeqChunks
    # TODO: Switch to Memory when possible
    # That also necessitates a first_index field
    v::Vector{UInt64}
    last_index::UInt
end

first_state(::LongSeqChunks) = UInt(1)

@inline function Base.iterate(it::LongSeqChunks, state::UInt=UInt(1))
    state > it.last_index && return nothing
    iter_inbounds(it, state)
end

@inline function iter_inbounds(it::LongSeqChunks, state::UInt)
    (@inbounds(it.v[state]), state + 1)
end

Base.eltype(::Type{LongSeqChunks}) = UInt64
Base.length(it::LongSeqChunks) = it.last_index % Int

function iter_chunks(seq::LongSequence)::Tuple{LongSeqChunks, Tuple{UInt64, UInt8}}
    v = seq.data
    seqlen = length(seq)
    n_bits = (seqlen * bits_per_symbol(Alphabet(seq))) % UInt
    remaining_bits = (n_bits % UInt(64)) % UInt8
    lst = index(lastbitindex(seq))
    len = !iszero(seqlen) * (lst - !iszero(remaining_bits))
    start = if iszero(remaining_bits)
        (UInt64(0), 0x00)
    else
        (@inbounds(v[lst]), remaining_bits)
    end
    (
        LongSeqChunks(v, len % UInt),
        start
    )
end

struct SubSeqChunks
    # TODO: Switch to Memory where possible
    v::Vector{UInt64}
    # First index to load from
    first_index::UInt
    # Last index to load from (except index + 1, if offset if not zero)
    last_index::UInt
    offset::UInt8
end

@inline function first_state(x::SubSeqChunks)
    x.first_index > x.last_index && return (x.last_index + 1, zero(UInt64))
    remaining = if iszero(x.offset)
        zero(UInt64)
    else
        @inbounds(x.v[x.first_index - 1]) >> (x.offset & 0x3f)
    end
    (x.first_index, remaining)
end

Base.eltype(::Type{SubSeqChunks}) = UInt64
Base.length(it::SubSeqChunks) = (it.last_index - it.first_index) % Int + 1

@inline function Base.iterate(it::SubSeqChunks, state::Tuple{UInt, UInt64}=first_state(it))
    first(state) > it.last_index && return nothing
    iter_inbounds(it, state)
end

@inline function iter_inbounds(it::SubSeqChunks, state::Tuple{UInt, UInt64})
    (index, remaining_bits) = state
    new_bits = @inbounds(it.v[index])
    result = (new_bits << ((0x40 - it.offset) & 0x3f)) | remaining_bits
    next_remaining = iszero(it.offset) ? zero(UInt64) : new_bits >> (it.offset & 0x3f)
    (result, (index + 1, next_remaining))
end


function iter_chunks(seq::LongSubSeq)::Tuple{SubSeqChunks, Tuple{UInt64, UInt8}}
    len = div(length(seq) % UInt, symbols_per_data_element(seq))
    bps = bits_per_symbol(Alphabet(seq))
    fbi = firstbitindex(seq)
    off = offset(fbi) % UInt8
    first_index = index(fbi) + !iszero(off)
    last_index = first_index + len - 1
    n_bits = (length(seq) * bps) % UInt
    remaining_bits = (n_bits % UInt(64)) % UInt8
    remaining_elems = div(remaining_bits, bps % UInt)
    v = seq.data
    start = if iszero(remaining_bits)
        (UInt64(0), 0x00)
    else
        lbi = lastbitindex(seq)
        bits_in_last = offset(lbi) + bps
        rsh = offset(bitindex(seq, lastindex(seq) - remaining_elems + 1))
        chunk = if bits_in_last < remaining_bits
            chunk = @inbounds(v[index(lbi) - 1]) >> (rsh & 63)
            chunk | @inbounds(v[index(lbi)]) << ((64 - rsh) & 63)
        else
            @inbounds(v[index(lbi)]) >> (rsh & 63)
        end
        (chunk, remaining_bits)
    end
    it = SubSeqChunks(v, first_index % UInt, last_index % UInt, off)
    (it, start)
end

struct PairedChunkIterator{A, B}
    a::A
    b::B
end

function iter_chunks(a::BioSequence, b::BioSequence)
    length(a) == length(b) || throw("Sequences must have same length")
    (ita, (cha, rema)) = iter_chunks(a)
    (itb, (chb, _)) = iter_chunks(b)
    (
        PairedChunkIterator{typeof(ita), typeof(itb)}(ita, itb),
        (cha, chb, rema),
    )
end

Base.length(it::PairedChunkIterator) = length(it.a)
Base.eltype(::Type{<:PairedChunkIterator}) = NTuple{2, UInt64}

@inline function Base.iterate(
    it::PairedChunkIterator,
    state=(first_state(it.a), first_state(it.b))
)
    a = iterate(it.a, first(state))
    isnothing(a) && return nothing
    b = iter_inbounds(it.b, last(state))
    ((first(a), first(b)), (last(a), last(b)))
end
