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
    (@inbounds(it.v[state]), state + one(state))
end

Base.eltype(::Type{LongSeqChunks}) = UInt64
Base.length(it::LongSeqChunks) = it.last_index % Int

function iter_chunks(seq::LongSequence)::Tuple{LongSeqChunks, Tuple{UInt64, UInt8}}
    v = seq.data
    n_bits = (length(seq) * bits_per_symbol(Alphabet(seq))) % UInt
    remaining_bits = (n_bits % UInt(64)) % UInt8
    lst = index(lastbitindex(seq))
    start = if iszero(remaining_bits)
        (UInt64(0), 0x00)
    else
        (@inbounds(v[lst]), remaining_bits)
    end
    (
        LongSeqChunks(v, (lst - !iszero(remaining_bits)) % UInt),
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

first_state(x::SubSeqChunks) = x.first_index

Base.eltype(::Type{SubSeqChunks}) = UInt64
Base.length(it::SubSeqChunks) = (it.last_index - it.first_index) % Int + 1

@inline function Base.iterate(it::SubSeqChunks, state::UInt=first_state(it))
    state > it.last_index && return nothing
    iter_inbounds(it, state)
end

@inline function iter_inbounds(it::SubSeqChunks, state::UInt)
    chunk = @inbounds(it.v[state]) >> (it.offset & 0x3f)
    chunk |= @inbounds(it.v[state + !iszero(it.offset)]) << ((0x40 - it.offset) & 0x3f)
    (chunk, state + one(state))
end

function iter_chunks(seq::LongSubSeq)::Tuple{SubSeqChunks, Tuple{UInt64, UInt8}}
    len = div(length(seq) % UInt, symbols_per_data_element(seq))
    bps = bits_per_symbol(Alphabet(seq))
    fbi = firstbitindex(seq)
    fst = index(fbi)
    off = offset(fbi) % UInt8
    lst = fst + len - 1
    n_bits = (length(seq) * bps) % UInt
    remaining_bits = (n_bits % UInt(64)) % UInt8
    remaining_elems = div(remaining_bits, bps % UInt)
    v = seq.data
    start = if iszero(remaining_bits)
        (UInt64(0), 0x00)
    else
        # TODO: This last cuhnk can also contain bits from the second-to-last index
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
    it = SubSeqChunks(v, fst % UInt, lst % UInt, off)
    (it, start)
end

struct PairedChunkIterator{A, B}
    a::A
    b::B
end

function iter_chunks(a::BioSequence, b::BioSequence)
    length(a) == length(b) || throw("Sequences must have same length")
    (ita, (cha, rema)) = iter_chunks(a)
    (itb, (chb, remb)) = iter_chunks(b)
    (
        PairedChunkIterator{typeof(ita), typeof(itb)}(ita, itb),
        (cha, chb, rema),
    )
end

Base.length(it::PairedChunkIterator) = length(it.a)
Base.eltype(::Type{<:PairedChunkIterator}) = UInt64

@inline function Base.iterate(it::PairedChunkIterator, state::UInt=first_state(it.a))
    a = iterate(it.a, state)
    isnothing(a) && return nothing
    b = iter_inbounds(it.b, state)
    ((first(a), first(b)), state + one(state))
end

@inline function counter_2seq(tail::Function, body::Function, a::SeqOrView, b::SeqOrView)
    (it, (ch1, ch2, rm)) = @inline iter_chunks(a, b)
    y = @inline tail(ch1, ch2, rm)
    for (a, b) in it
        y += @inline body(a, b)
    end
    y 
end

@inline function counter_1seq(tail::Function, body::Function, seq::SeqOrView)
    (it, (chunk, rm)) = @inline iter_chunks(seq)
    y = @inline tail(chunk, rm)
    for i in it
        y += @inline body(i)
    end
    y
end

# TODO: We should change this API in a future biojulia release:
# * This does not conform to the count interface, so delete the count overloads
# * Remove most two-seq methods, these are nonsensical

"""
    gc_content(seq::BioSequence) -> Float64

Calculate GC content of `seq`, i.e. the number of symbols that is `DNA_C`, `DNA_G`,
`DNA_C` or `DNA_G` divided by the length of the sequence.

# Examples
```jldoctest
julia> gc_content(dna"AGCTA")
0.4

julia> gc_content(rna"UAGCGA")
0.5
```
"""
gc_content(seq::NucleotideSeq) = count(isGC, seq) / length(seq)

# Aliases
"""
    mismatches(a::BioSequence, b::BioSequences) -> Int

Count the number of positions in where `a` and `b` are equal.
If `b` is given, and the length of `a` and `b` differ, look only at the indices
of the shorter sequence.

# Examples
```jldoctest
julia> matches(dna"TAGCTA", dna"TACCTA")
5

julia> matches(dna"AACA", dna"AAG")
2
```
"""
mismatches(a::BioSequence, b::BioSequence) = count(!=, a, b)

"""
    mismatches(a::BioSequence, b::BioSequences) -> Int

Count the number of positions in where `a` and `b` differ.
If `b` is given, and the length of `a` and `b` differ, look only at the indices
of the shorter sequence.

# Examples
```jldoctest
julia> mismatches(dna"TAGCTA", dna"TACCTA")
1

julia> mismatches(dna"AACA", dna"AAG")
1
```
"""
matches(a::BioSequence, b::BioSequence) = count(==, a, b)

"""
    n_gaps(a::BioSequence, [b::BioSequence]) -> Int

Count the number of positions where `a` (or `b`, if present) have gaps.
If `b` is given, and the length of `a` and `b` differ, look only at the indices
of the shorter sequence.

# Examples
```jldoctest
julia> n_gaps(dna"--TAC-WN-ACY")
4

julia> n_gaps(dna"TC-AC-", dna"-CACG")
2
```
"""
function n_gaps end

n_gaps(seq::BioSequence) = count(isgap, seq)
n_gaps(a::BioSequence, b::BioSequence) = count(isgap, a, b)

"""
    n_ambiguous(a::BioSequence, [b::BioSequence]) -> Int

Count the number of positions where `a` (or `b`, if present) have ambigious symbols.
If `b` is given, and the length of `a` and `b` differ, look only at the indices
of the shorter sequence.

# Examples
```jldoctest
julia> n_ambiguous(dna"--TAC-WN-ACY")
3

julia> n_ambiguous(rna"UAYWW", rna"UAW")
1
```
"""
function n_ambiguous end

n_ambiguous(seq::BioSequence) = count(isambiguous, seq)
n_ambiguous(a::BioSequence, b::BioSequence) = count(isambiguous, a, b)

"""
    n_certain(a::BioSequence, [b::BioSequence]) -> Int

Count the number of positions where `a` (and `b`, if present) have certain (i.e. non-ambigous
and non-gap) symbols.
If `b` is given, and the length of `a` and `b` differ, look only at the indices
of the shorter sequence.

# Examples
```jldoctest
julia> n_certain(dna"--TAC-WN-ACY")
5

julia> n_certain(rna"UAYWW", rna"UAW")
2
```
"""
function n_certain end

n_certain(seq::BioSequence) = count(iscertain, seq)
n_certain(a::BioSequence, b::BioSequence) = count(iscertain, a, b)

# Fallback implementation
# (Urgh, we really shouldn't do these kinds of overloads of Base)
function count_naive(pred, seqa::BioSequence, seqb::BioSequence)
    n = 0
    for (i, j) in zip(seqa, seqb)
        n += pred(i, j)::Bool
    end
    return n
end

Base.count(f::Function, a::BioSequence, b::BioSequence) = count_naive(f, a, b)

# Overloads of count
const FourBit = SeqOrView{<:NucleicAcidAlphabet{4}}
const TwoBit = SeqOrView{<:NucleicAcidAlphabet{2}}

# Trivial overloads
Base.count(::typeof(iscertain), s::BioSequence{<:NucleicAcidAlphabet{2}}) = length(s)
Base.count(::typeof(isambiguous), s::BioSequence{<:NucleicAcidAlphabet{2}}) = 0
Base.count(::typeof(isgap), s::BioSequence{<:NucleicAcidAlphabet{2}}) = 0

function Base.count(
    ::typeof(iscertain),
    a::BioSequence{<:NucleicAcidAlphabet{2}},
    b::BioSequence{<:NucleicAcidAlphabet{2}},
)
    min(length(a), length(b))
end

function Base.count(
    ::typeof(isgap),
    a::BioSequence{<:NucleicAcidAlphabet{2}},
    b::BioSequence{<:NucleicAcidAlphabet{2}},
)
    0
end

function Base.count(
    ::typeof(isambiguous),
    a::BioSequence{<:NucleicAcidAlphabet{2}},
    b::BioSequence{<:NucleicAcidAlphabet{2}},
)
    0
end

# TODO: For ambiguous, certain, isgap with mixed 2/4 bit nucl, we can take shortcuts
# TODO: Also implement this for AA sequences

function Base.count(::typeof(==), a::BioSequence, b::BioSequence)
    min(length(a), length(b)) - count(!=, a, b)
end

# Other overloads
function Base.count(
    ::typeof(!=),
    a::TwoBit,
    b::TwoBit,
)
    tail = (ch1, ch2, rm) -> begin
        mask = UInt64(1) << (rm & 63) - 1
        count_nonzero_bitpairs((ch1 & mask) ⊻ (ch2 & mask))
    end
    body = (i, j) -> count_nonzero_bitpairs(i ⊻ j)
    counter_2seq(tail, body, a, b)
end

function Base.count(
    ::typeof(!=),
    a::FourBit,
    b::FourBit,
)
    tail = (ch1, ch2, rm) -> begin
        mask = UInt64(1) << (rm & 63) - 1
        count_nonzero_nibbles((ch1 & mask) ⊻ (ch2 & mask))
    end
    body = (i, j) -> count_nonzero_nibbles(i ⊻ j)
    counter_2seq(tail, body, a, b)
end

function Base.count(::typeof(isGC), seq::Union{TwoBit, FourBit})
    tail = (chunk, rm) -> begin
        mask = UInt64(1) << (rm & 63) - 1
        gc_bitcount(chunk & mask, Alphabet(seq))
    end
    body = i -> gc_bitcount(i, Alphabet(seq))
    counter_1seq(tail, body, seq)
end

function Base.count(::typeof(isgap), seq::FourBit)
    tail = (chunk, rm) -> begin
        mask = ~(UInt64(1) << (rm & 63) - 1)
        count_0000_nibbles(chunk | mask)
    end
    counter_1seq(tail, count_0000_nibbles, seq)
end

function Base.count(::typeof(isgap), a::FourBit, b::FourBit)
    tail = (ch1, ch2, rm) -> begin
        mask = ~(UInt64(1) << (rm & 63) - 1)
        gap_bitcount(ch1 | mask, ch2 | mask, Alphabet(a))
    end
    body = (i, j) -> gap_bitcount(i, j, Alphabet(a))
    counter_2seq(tail, body, a, b)
end

function Base.count(::typeof(isambiguous), seq::FourBit)
    tail = (chunk, rm) -> begin
        mask = (UInt64(1) << (rm & 63) - 1)
        ambiguous_bitcount(chunk & mask, Alphabet(seq))
    end
    body = i -> ambiguous_bitcount(i, Alphabet(seq))
    counter_1seq(tail, body, seq)
end

function Base.count(::typeof(isambiguous), a::FourBit, b::FourBit)
    tail = (ch1, ch2, rm) -> begin
        mask = (UInt64(1) << (rm & 63) - 1)
        ambiguous_bitcount(ch1 & mask, ch2 & mask, Alphabet(seq))
    end
    body = (i, j) -> ambiguous_bitcount(i, j, Alphabet(seq))
    counter_2seq(tail, body, a, b)
end

function Base.count(::typeof(iscertain), seq::FourBit)
    tail = (chunk, rm) -> begin
        mask = (UInt64(1) << (rm & 63) - 1)
        certain_bitcount(chunk & mask, Alphabet(seq))
    end
    body = i -> certain_bitcount(i, Alphabet(seq))
    counter_1seq(tail, body, seq)
end

function Base.count(::typeof(isambiguous), a::FourBit, b::FourBit)
    tail = (ch1, ch2, rm) -> begin
        mask = (UInt64(1) << (rm & 63) - 1)
        certain_bitcount(ch1 & mask, ch2 & mask, Alphabet(seq))
    end
    body = (i, j) -> certain_bitcount(i, j, Alphabet(seq))
    counter_2seq(tail, body, a, b)
end
 
