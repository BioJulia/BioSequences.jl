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

# This returns (head, body, tail), where:
# - head and tail are Tuple{UInt64, UInt8}, with a coding element and the number
#   of coding bits in that element. Head is the partial coding element before any
#   full elements, and tail is the partial after any coding elements.
#   If head or tail is empty, the UInt8 is set to zero. By definition, it can be
#   at most set to 63.
#   If the sequence is composed of only one partial element, tail is nonempty
#   and head is empty.
# - body is a Tuple{UInt, UInt} with the (start, stop) indices of coding elements.
#   If stop < start, there are no such elements.
# TODO: The body should probably be a MemoryView in 1.11
function parts(seq::LongSequence)
    # LongSequence never has coding bits before the first chunks
    head = (zero(UInt64), zero(UInt8))
    len = length(seq)
    # Shortcut to prevent annoying edge cases in the rest of the code
    if iszero(len)
        return (head, (UInt(1), UInt(0)), (zero(UInt64), zero(UInt8)))
    end
    lastbitindex(seq)
    bits_in_tail = (offset(bitindex(seq, len + 1)) % UInt8) & 0x3f
    lbi = bitindex(seq, len)
    lbii = index(lbi)
    tail = if iszero(bits_in_tail)
        head
    else
        (@inbounds(seq.data[lbii]), bits_in_tail)
    end
    # If we have bits in the tail, then clearly those bits means the last bitindex
    # points to one past the last full chunk
    body = (UInt(1), (lbii - !iszero(bits_in_tail)) % UInt)
    (head, body, tail)
end

function parts(seq::LongSubSeq)
    data = seq.data
    zero_end = (zero(UInt64), zero(UInt8))
    len = length(seq)
    # Again: Avoid annoying edge cases later
    if iszero(len)
        return (zero_end, (UInt(1), UInt(0)), zero_end)
    end
    lastbitindex(seq)
    lbi = bitindex(seq, len)
    lbii = index(lbi)
    fbi = firstbitindex(seq)
    fbii = index(fbi)
    bits_in_head_naive = (((64 - offset(fbi)) % UInt8) & 0x3f)
    # If first and last chunk index is the same, there are actually zero
    # bits in head, as they are all in the tail
    bits_in_head = bits_in_head_naive * (lbii != fbii)
    # For the head, there are some uncoding lower bits. We need to shift
    # the head right with this number.
    head_shift = ((0x40 - bits_in_head_naive) & 0x3f)
    head = if iszero(bits_in_head)
        zero_end
    else
        chunk = @inbounds(data[fbii]) >> head_shift
        (chunk, bits_in_head)
    end
    # However, if last and first chunk index is the same, there is no head
    # chunk, and thus no head chunk to shift, but the TAIL chunk may not have coding bits at the lowest
    # position.
    tail_shift = (head_shift * (lbii == fbii)) & 63
    bits_in_tail = (offset(bitindex(seq, len + 1)) % UInt8) & 0x3f
    bits_in_tail -= tail_shift % UInt8
    tail = if iszero(bits_in_tail)
        zero_end
    else
        (@inbounds(data[lbii]) >> tail_shift, bits_in_tail)
    end
    body = (fbii + !iszero(bits_in_head), lbii - !iszero(bits_in_tail))
    (head, body, tail)
end