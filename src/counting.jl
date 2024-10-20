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

first_state(x::SubSeqChunks) = x.first_index

Base.eltype(::Type{SubSeqChunks}) = UInt64
Base.length(it::SubSeqChunks) = (it.last_index - it.first_index) % Int + 1

@inline function Base.iterate(it::SubSeqChunks, state::UInt=first_state(it))
    state > it.last_index && return nothing
    iter_inbounds(it, state)
end

# TODO: We could have one memory load per iteration, and store unused bits in the state
# for the next iteration
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
    state::Tuple{UInt, UInt}=(first_state(it.a), first_state(it.b))
)
    a = iterate(it.a, first(state))
    isnothing(a) && return nothing
    b = iter_inbounds(it.b, last(state))
    ((first(a), first(b)), (first(state) + 1, last(state) + 1))
end

@inline function counter_2seq(tail::Function, body::Function, a::SeqOrView, b::SeqOrView)
    # TODO: Remove these first lines at the API change
    minlength = min(length(a), length(b))
    (a, b) = (view(a, 1:minlength), view(b, 1:minlength))
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
gc_content(seq::NucleotideSeq) = n_gc(seq) / length(seq)
Base.count(::typeof(isGC), seq::NucleotideSeq) = n_gc(seq)

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
Base.count(::typeof(!=), a::BioSequence, b::BioSequence) = mismatches(a, b)

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
Base.count(::typeof(==), a::BioSequence, b::BioSequence) = matches(a, b)

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

Base.count(::typeof(isgap), seq::BioSequence) = n_gaps(seq)
Base.count(::typeof(isgap), a::BioSequence, b::BioSequence) = n_gaps(a, b)

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

Base.count(::typeof(isambiguous), seq::BioSequence) = n_ambiguous(seq)
Base.count(::typeof(isambiguous), a::BioSequence, b::BioSequence) = n_ambiguous(a, b)

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

Base.count(::typeof(iscertain), seq::BioSequence) = n_certain(seq)
Base.count(::typeof(iscertain), a::BioSequence, b::BioSequence) = n_certain(a, b)

# Fallback implementation
# (Urgh, we really shouldn't do these kinds of overloads of Base)
function count_naive(pred, seqa::BioSequence, seqb::BioSequence)
    n = 0
    for (i, j) in zip(seqa, seqb)
        n += pred(i, j)::Bool
    end
    return n
end

function count_naive(pred, seq::BioSequence)
    n = 0
    for i in seq
        n += pred(i)::Bool
    end
    return n
end

Base.count(f::Function, a::BioSequence, b::BioSequence) = count_naive(f, a, b)

mismatches(a::BioSequence, b::BioSequence) = count_naive(!=, a, b)

n_ambiguous(seq::BioSequence) = count_naive(isambiguous, seq)
function n_ambiguous(a::BioSequence, b::BioSequence)
    count_naive((x, y) -> isambiguous(x) | isambiguous(y), a, b)
end

n_gaps(seq::BioSequence) = count_naive(isgap, seq)
function n_gaps(a::BioSequence, b::BioSequence)
    count_naive((x, y) -> isgap(x) | isgap(y), a, b)
end

n_certain(seq::BioSequence) = count_naive(iscertain, seq)
function n_certain(a::BioSequence, b::BioSequence)
    count_naive((x, y) -> iscertain(x) & iscertain(y), a, b)
end
n_gc(seq::BioSequence) = count_naive(isGC, seq)

const FourBit = SeqOrView{<:NucleicAcidAlphabet{4}}
const TwoBit = SeqOrView{<:NucleicAcidAlphabet{2}}

const TwoBitSeq = BioSequence{<:NucleicAcidAlphabet{2}}
const FourBitSeq = BioSequence{<:NucleicAcidAlphabet{4}}

# Trivial overloads
n_certain(s::TwoBitSeq) = length(s)
n_ambiguous(s::TwoBitSeq) = 0
n_gaps(s::TwoBitSeq) = 0

n_certain(a::TwoBitSeq, b::TwoBitSeq) = min(length(a), length(b))
n_gaps(::TwoBitSeq, ::TwoBitSeq) = 0
n_ambiguous(::TwoBitSeq, ::TwoBitSeq) = 0

function matches(a::BioSequence, b::BioSequence)
    min(length(a), length(b)) - mismatches(a, b)
end

# Other overloads

# Mismatches is special because DNA_T != RNA_U despite them having the same encoding,
# so we need to use the fallback
# TODO: Or have a serious hack counting 1000 nibbles / 11 bitpairs.
function mismatches(a::SeqOrView{A}, b::SeqOrView{A}) where {A <: NucleicAcidAlphabet{2}}
    tail = (ch1, ch2, rm) -> begin
        mask = UInt64(1) << (rm & 63) - 1
        count_nonzero_bitpairs((ch1 & mask) ⊻ (ch2 & mask))
    end
    body = (i, j) -> count_nonzero_bitpairs(i ⊻ j)
    counter_2seq(tail, body, a, b)
end

function mismatches(a::SeqOrView{A}, b::SeqOrView{A}) where {A <: NucleicAcidAlphabet{4}}
    tail = (ch1, ch2, rm) -> begin
        mask = UInt64(1) << (rm & 63) - 1
        count_nonzero_nibbles((ch1 & mask) ⊻ (ch2 & mask))
    end
    body = (i, j) -> count_nonzero_nibbles(i ⊻ j)
    counter_2seq(tail, body, a, b)
end

function n_gc(seq::Union{TwoBit, FourBit})
    tail = (chunk, rm) -> begin
        mask = UInt64(1) << (rm & 63) - 1
        gc_bitcount(chunk & mask, Alphabet(seq))
    end
    body = i -> gc_bitcount(i, Alphabet(seq))
    counter_1seq(tail, body, seq)
end

function n_gaps(seq::FourBit)
    tail = (chunk, rm) -> begin
        mask = ~(UInt64(1) << (rm & 63) - 1)
        count_0000_nibbles(chunk | mask)
    end
    counter_1seq(tail, count_0000_nibbles, seq)
end

function n_gaps(a::FourBit, b::FourBit)
    tail = (ch1, ch2, rm) -> begin
        mask = ~(UInt64(1) << (rm & 63) - 1)
        gap_bitcount(ch1 | mask, ch2 | mask, Alphabet(a))
    end
    body = (i, j) -> gap_bitcount(i, j, Alphabet(a))
    counter_2seq(tail, body, a, b)
end

function n_ambiguous(seq::FourBit)
    tail = (chunk, rm) -> begin
        mask = (UInt64(1) << (rm & 63) - 1)
        ambiguous_bitcount(chunk & mask, Alphabet(seq))
    end
    body = i -> ambiguous_bitcount(i, Alphabet(seq))
    counter_1seq(tail, body, seq)
end

function n_ambiguous(a::FourBit, b::FourBit)
    tail = (ch1, ch2, rm) -> begin
        mask = (UInt64(1) << (rm & 63) - 1)
        ambiguous_bitcount(ch1 & mask, ch2 & mask, Alphabet(a))
    end
    body = (i, j) -> ambiguous_bitcount(i, j, Alphabet(a))
    counter_2seq(tail, body, a, b)
end

function n_certain(seq::FourBit)
    tail = (chunk, rm) -> begin
        mask = (UInt64(1) << (rm & 63) - 1)
        certain_bitcount(chunk & mask, Alphabet(seq))
    end
    body = i -> certain_bitcount(i, Alphabet(seq))
    counter_1seq(tail, body, seq)
end

function n_certain(a::FourBit, b::FourBit)
    tail = (ch1, ch2, rm) -> begin
        mask = (UInt64(1) << (rm & 63) - 1)
        certain_bitcount(ch1 & mask, ch2 & mask, Alphabet(a))
    end
    body = (i, j) -> certain_bitcount(i, j, Alphabet(a))
    counter_2seq(tail, body, a, b)
end

n_ambiguous(a::TwoBit, b::FourBit) = n_ambiguous(b, a)
n_certain(a::TwoBit, b::FourBit) = n_certain(b, a)
n_gaps(a::TwoBit, b::FourBit) = n_gaps(b, a)

function n_ambiguous(a::FourBit, b::TwoBit)
    n_ambiguous(view(a, 1:min(length(a), length(b))))
end

function n_certain(a::FourBit, b::TwoBit)
    minlength = min(length(a), length(b))
    n_certain(view(a, 1:minlength))
end

function n_gaps(a::FourBit, b::TwoBit)
    n_gaps(view(a, 1:min(length(a), length(b))))
end

# TODO: Also implement this for AA sequences
# mismatches
# n_gaps
# n_ambiguous
# n_certain


# Gap is 0x1b
# Ambig: 0x16:0x19
# Stop: 0x1a