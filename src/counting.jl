trunc_seq(x::LongSequence, len::Int) = typeof(x)(x.data, len % UInt)
trunc_seq(x::LongSubSeq, len::Int) = typeof(x)(x.data, first(x.part):first(x.part)+len-1)

@inline function counter_2seq(tail::Function, body::Function, a::SeqOrView, b::SeqOrView)
    # TODO: Remove these first lines at the API change
    minlength = min(length(a), length(b))
    (a, b) = (trunc_seq(a, minlength), trunc_seq(b, minlength))
    (it, (ch1, ch2, rm)) = @inline iter_chunks(a, b)
    y = @inline tail(ch1, ch2, rm)
    for (a, b) in it
        y += @inline body(a, b)
    end
    y 
end

@inline function counter_1seq(partial::Function, body::Function, seq::SeqOrView)
    ((head, head_bits), (chunk_start, chunk_end), (tail, tail_bits)) = parts(seq)
    y = @inline partial(head, head_bits) + @inline partial(tail, tail_bits)
    data = seq.data
    @inbounds for i in chunk_start:chunk_end
        y += @inline body(data[i])
    end
    y
    y
end

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
gc_content(seq::NucleotideSeq) = _n_gc(seq) / length(seq)
Base.count(::typeof(isGC), seq::NucleotideSeq) = _n_gc(seq)

# Aliases
"""
    matches(a::BioSequence, b::BioSequences) -> Int

Count the number of positions in where `a` and `b` are equal.
If `b` is given, and the length of `a` and `b` differ, look only at the indices
of the shorter sequence.
This function does not provide any special handling of ambiguous symbols,
so e.g. `DNA_A` does not match `DNA_N`.

!!! warning
    Passing in two sequences with differing lengths is deprecated. In a future,
    breaking release of BioSequences, this will error.

# Examples
```jldoctest
julia> matches(dna"TAWNNA", dna"TACCTA")
3

julia> matches(dna"AACA", dna"AAG")
2
```
"""
function matches end

@deprecate(
    Base.count(::typeof(!=), a::BioSequence, b::BioSequence),
    count(splat(!=), zip(a, b)),
    false
)

"""
    mismatches(a::BioSequence, b::BioSequences) -> Int

Count the number of positions in where `a` and `b` differ.
If `b` is given, and the length of `a` and `b` differ, look only at the indices
of the shorter sequence.
This function does not provide any special handling of ambiguous symbols,
so e.g. `DNA_A` does not match `DNA_N`.

!!! warning
    Passing in two sequences with differing lengths is deprecated. In a future,
    breaking release of BioSequences, this will error.

# Examples
```jldoctest
julia> mismatches(dna"TAGCTA", dna"TACNTA")
2

julia> mismatches(dna"AACA", dna"AAG")
1
```
"""
function mismatches end

@deprecate(
    Base.count(::typeof(==), a::BioSequence, b::BioSequence),
    count(splat(==), zip(a, b)),
    false
)

"""
    n_gaps(a::BioSequence, [b::BioSequence]) -> Int

Count the number of positions where `a` (or `b`, if present) have gaps.
If `b` is given, and the length of `a` and `b` differ, look only at the indices
of the shorter sequence.

!!! warning
    Passing in two sequences is deprecated. In a future, breaking release of
    BioSequences, this will throw a `MethodError`

# Examples
```jldoctest
julia> n_gaps(dna"--TAC-WN-ACY")
4

julia> n_gaps(dna"TC-AC-", dna"-CACG")
2
```
"""
function n_gaps end

Base.count(::typeof(isgap), seq::BioSequence) = _n_gaps(seq)

@deprecate(
    Base.count(::typeof(isgap), a::BioSequence, b::BioSequence),
    count(((i,j),) -> isgap(i) | isgap(j), zip(a, b)),
    false
)

"""
    n_ambiguous(a::BioSequence, [b::BioSequence]) -> Int

Count the number of positions where `a` (or `b`, if present) have ambigious symbols.
If `b` is given, and the length of `a` and `b` differ, look only at the indices
of the shorter sequence.
Gaps are not ambigous.

!!! warning
    Passing in two sequences is deprecated. In a future, breaking release of
    BioSequences, this will throw a `MethodError`


# Examples
```jldoctest
julia> n_ambiguous(dna"--TAC-WN-ACY")
3

julia> n_ambiguous(rna"UAYWW", rna"UAW")
1
```
"""
function n_ambiguous end

Base.count(::typeof(isambiguous), seq::BioSequence) = _n_ambiguous(seq)

@deprecate(
    Base.count(::typeof(isambiguous), a::BioSequence, b::BioSequence),
    count(((i,j),) -> isambiguous(i) | isambiguous(j), zip(a, b)),
    false
)

"""
    n_certain(a::BioSequence, [b::BioSequence]) -> Int

Count the number of positions where `a` (and `b`, if present) have certain (i.e. non-ambigous
and non-gap) symbols.
If `b` is given, and the length of `a` and `b` differ, look only at the indices
of the shorter sequence.
Gaps are not certain.

!!! warning
    Passing in two sequences is deprecated. In a future, breaking release of
    BioSequences, this will throw a `MethodError`


# Examples
```jldoctest
julia> n_certain(dna"--TAC-WN-ACY")
5

julia> n_certain(rna"UAYWW", rna"UAW")
2
```
"""
function n_certain end

Base.count(::typeof(iscertain), seq::BioSequence) = _n_certain(seq)

@deprecate(
    Base.count(::typeof(iscertain), a::BioSequence, b::BioSequence),
    count(((i,j),) -> iscertain(i) & iscertain(j), zip(a, b)),
    false
)

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

@deprecate(
    Base.count(f::Function, a::BioSequence, b::BioSequence),
    count(splat(f), a, b),
    false
)

@noinline function depwarn_lengths()
    Base.depwarn(
        "Calling `mismatches` with sequences of differing lengths is deprecated,
        and will throw an exception in a future release of BioSequences.",
        :mismatches
    )
end

function mismatches_inner(a::BioSequence, b::BioSequence)
    length(a) == length(b) || depwarn_lengths()
    count_naive(!=, a, b)
end

# Compile-time check with conservative fallback to `false`
# for unknown alphabets
Base.@assume_effects :foldable function known_disjoint(a::KNOWN_ALPHABETS, b::KNOWN_ALPHABETS)
    for i in a, j in b
        i == j && return false
    end
    true
end

known_disjoint(::Alphabet, ::Alphabet) = false

# For the known alphabets (vast majojrity of cases), we can tell statically
# if the alphabets don't overlap, e.g. DNAAlphabet{2} and RNAAlphabet{4}.
# In this case, the result will always be zero.
function mismatches(a::BioSequence{A}, b::BioSequence{B}) where {A, B}
    known_disjoint(A(), B()) ? min(length(a), length(b)) : mismatches_inner(a, b)
end

_n_ambiguous(seq::BioSequence) = count_naive(isambiguous, seq)

_n_certain(seq::BioSequence) = count_naive(iscertain, seq)

_n_gc(seq::BioSequence) = count_naive(isGC, seq)

# Trivial overloads
_n_certain(s::TwoBitSeq) = length(s)
_n_ambiguous(s::TwoBitSeq) = 0
_n_gaps(s::TwoBitSeq) = 0

function matches(a::BioSequence, b::BioSequence)
    min(length(a), length(b)) - mismatches(a, b)
end

_n_gaps(s::BioSequence) = count_symbol(s, gap(eltype(Alphabet(s))))

# Other overloads
function mismatches_inner(a::SeqOrView{A}, b::SeqOrView{A}) where {A <: NucleicAcidAlphabet{2}}
    length(a) == length(b) || depwarn_lengths()
    tail = (ch1, ch2, rm) -> begin
        mask = UInt64(1) << (rm & 63) - 1
        count_nonzero_bitpairs((ch1 & mask) ⊻ (ch2 & mask))
    end
    body = (i, j) -> count_nonzero_bitpairs(i ⊻ j)
    counter_2seq(tail, body, a, b)
end

function mismatches_inner(a::SeqOrView{A}, b::SeqOrView{A}) where {A <: NucleicAcidAlphabet{4}}
    length(a) == length(b) || depwarn_lengths()
    tail = (ch1, ch2, rm) -> begin
        mask = UInt64(1) << (rm & 63) - 1
        count_nonzero_nibbles((ch1 & mask) ⊻ (ch2 & mask))
    end
    body = (i, j) -> count_nonzero_nibbles(i ⊻ j)
    counter_2seq(tail, body, a, b)
end

function _n_gc(seq::Union{TwoBit, FourBit})
    tail = (chunk, rm) -> begin
        mask = UInt64(1) << (rm & 63) - 1
        gc_bitcount(chunk & mask, Alphabet(seq))
    end
    body = i -> gc_bitcount(i, Alphabet(seq))
    counter_1seq(tail, body, seq)
end

function _n_ambiguous(seq::FourBit)
    tail = (chunk, rm) -> begin
        mask = (UInt64(1) << (rm & 63) - 1)
        ambiguous_bitcount(chunk & mask, Alphabet(seq))
    end
    body = i -> ambiguous_bitcount(i, Alphabet(seq))
    counter_1seq(tail, body, seq)
end

function _n_certain(seq::FourBit)
    tail = (chunk, rm) -> begin
        mask = (UInt64(1) << (rm & 63) - 1)
        certain_bitcount(chunk & mask, Alphabet(seq))
    end
    body = i -> certain_bitcount(i, Alphabet(seq))
    counter_1seq(tail, body, seq)
end

function mismatches_inner(a::SeqOrView{AminoAcidAlphabet}, b::SeqOrView{AminoAcidAlphabet})
    length(a) == length(b) || depwarn_lengths()
    tail = (ch1, ch2, rm) -> begin
        mask = UInt64(1) << (rm & 63) - 1
        count_compared_bytes(!iszero, (ch1 & mask) ⊻ (ch2 & mask))
    end
    body = (i, j) -> count_compared_bytes(!iszero, i ⊻ j)
    counter_2seq(tail, body, a, b)
end

function _n_ambiguous(seq::SeqOrView{AminoAcidAlphabet})
    tail = (chunk, rm) -> begin
        mask = UInt64(1) << (rm & 63) - 1
        ambiguous_bitcount(chunk & mask, AminoAcidAlphabet())
    end
    body = i -> ambiguous_bitcount(i, AminoAcidAlphabet())
    counter_1seq(tail, body, seq)
end

function _n_certain(seq::SeqOrView{AminoAcidAlphabet})
    tail = (chunk, rm) -> begin
        mask = ~(UInt64(1) << (rm & 63) - 1)
        count_compared_bytes(<(0x16), chunk | mask)
    end
    body = i -> count_compared_bytes(<(0x16), i)
    counter_1seq(tail, body, seq)
end

function count_symbol(seq::BioSequence, sym::BioSymbol)
    # We encode here in order to throw an exception
    encode(Alphabet(seq), sym)
    n = 0
    for i in seq
        n += (i == sym)::Bool
    end
    n
end

function Base.count(
    pred::Base.Fix2{<:Union{typeof(==), typeof(isequal)}, <: BioSymbol},
    s::BioSequence
)
    count_symbol(s, pred.x)
end

function count_symbol(seq::FourBit, s::Union{RNA, DNA})
    pattern = encode(Alphabet(seq), s)::UInt64 * 0x1111111111111111
    tail = (chunk, rm) -> begin
        mask = UInt64(1) << (rm & 63) - 1
        masked = iszero(pattern) ? chunk | ~mask : mask & chunk
        count_0000_nibbles(masked ⊻ pattern)
    end
    body = i -> count_0000_nibbles(i ⊻ pattern)
    counter_1seq(tail, body, seq)
end

function count_symbol(seq::TwoBit, s::Union{RNA, DNA})
    pattern = encode(Alphabet(seq), s)::UInt64 * 0x5555555555555555
    tail = (chunk, rm) -> begin
        mask = UInt64(1) << (rm & 63) - 1
        masked = iszero(pattern) ? chunk | ~mask : mask & chunk
        count_00_bitpairs(masked ⊻ pattern)
    end
    body = i -> count_00_bitpairs(i ⊻ pattern)
    counter_1seq(tail, body, seq)
end

function count_symbol(seq::SeqOrView{AminoAcidAlphabet}, s::AminoAcid)
    byte = encode(Alphabet(seq), s) % UInt8
    tail = (chunk, rm) -> begin
        mask = UInt64(1) << (rm & 63) - 1
        masked = iszero(byte) ? chunk | ~mask : mask & chunk
        count_compared_bytes(==(byte), masked)
    end
    body = i -> count_compared_bytes(==(byte), i)
    counter_1seq(tail, body, seq)
end

## Deprecate weird two-arg methods
for (sym, op, joiner) in [
    (:n_ambiguous, :isambiguous, :|),
    (:n_gaps, :isgap, :|),
    (:n_certain, :iscertain, :&),
]
    for (A, B) in [
        (:BioSequence, :BioSequence),
        (:FourBit, :FourBit),
        (:FourBit, :TwoBit),
        (:TwoBit, :FourBit),
        (:TwoBit, :TwoBit),
    ]
        @eval(@deprecate(
            $(sym)(a::$(A), b::$(B)),
            count(((i,j),) -> $(joiner)($(op)(i), $(op)(j)), zip(a, b)),
        ))
    end
    for T in [:BioSequence, TwoBit, FourBit, :(SeqOrView{<:AminoAcidAlphabet})]
        @eval(@deprecate($(sym)(a::$(T)), count($(op), a)))
    end
end

