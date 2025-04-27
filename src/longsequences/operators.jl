
# Sequences to Matrix
# -------------------

"""
    seqmatrix(vseq::AbstractVector{BioSequence{A}}, major::Symbol) where {A<:Alphabet}

Construct a matrix of nucleotides or amino acids from a vector of `BioSequence`s.

If parameter `major` is set to `:site`, the matrix is created such that one
nucleotide from each sequence is placed in each column i.e. the matrix is laid
out in site-major order.
This means that iteration over one position of many sequences is efficient,
as julia arrays are laid out in column major order.

If the parameter `major` is set to `:seq`, the matrix is created such that each
sequence is placed in one column i.e. the matrix is laid out in sequence-major
order.
This means that iteration across each sequence in turn is efficient,
as julia arrays are laid out in column major order.

# Examples
```julia
julia> seqs = [dna"AAA", dna"TTT", dna"CCC", dna"GGG"]
4-element Array{BioSequences.BioSequence{BioSequences.DNAAlphabet{4}},1}:
 3nt DNA Sequence:
AAA
 3nt DNA Sequence:
TTT
 3nt DNA Sequence:
CCC
 3nt DNA Sequence:
GGG

julia> seqmatrix(seqs, :site)
4x3 Array{BioSequences.DNA,2}:
 DNA_A  DNA_A  DNA_A
 DNA_T  DNA_T  DNA_T
 DNA_C  DNA_C  DNA_C
 DNA_G  DNA_G  DNA_G

 julia> seqmatrix(seqs, :seq)
 3x4 Array{BioSequences.DNA,2}:
  DNA_A  DNA_T  DNA_C  DNA_G
  DNA_A  DNA_T  DNA_C  DNA_G
  DNA_A  DNA_T  DNA_C  DNA_G
```
"""
function seqmatrix(vseq::AbstractVector{LongSequence{A}}, major::Symbol) where {A<:Alphabet}
    nseqs = length(vseq)
    @assert nseqs > 0 throw(ArgumentError("Vector of BioSequence{$A} is empty."))
    nsites = length(vseq[1])
    @inbounds for i in 2:nseqs
        length(vseq[i]) == nsites || throw(ArgumentError("Sequences in vseq must be of same length"))
    end
    if major == :site
        mat = Matrix{eltype(A)}(undef, (nseqs, nsites))
        @inbounds for seq in 1:nseqs, site in 1:nsites
            mat[seq, site] = vseq[seq][site]
        end
        return mat
    elseif major == :seq
        mat = Matrix{eltype(A)}(undef, (nsites, nseqs))
        @inbounds for seq in 1:nseqs, site in 1:nsites
            mat[site, seq] = vseq[seq][site]
        end
        return mat
    else
        throw(ArgumentError("major must be :site or :seq"))
    end
end

"""
    seqmatrix(::Type{T}, vseq::AbstractVector{BioSequence{A}}, major::Symbol) where {T,A<:Alphabet}

Construct a matrix of `T` from a vector of `BioSequence`s.

If parameter `major` is set to `:site`, the matrix is created such that one
nucleotide from each sequence is placed in each column i.e. the matrix is laid
out in site-major order.
This means that iteration over one position of many sequences is efficient,
as julia arrays are laid out in column major order.

If the parameter `major` is set to `:seq`, the matrix is created such that each
sequence is placed in one column i.e. the matrix is laid out in sequence-major
order.
This means that iteration across each sequence in turn is efficient,
as julia arrays are laid out in column major order.

# Examples
```julia
julia> seqs = [dna"AAA", dna"TTT", dna"CCC", dna"GGG"]
4-element Array{BioSequences.BioSequence{BioSequences.DNAAlphabet{4}},1}:
 3nt DNA Sequence:
AAA
 3nt DNA Sequence:
TTT
 3nt DNA Sequence:
CCC
 3nt DNA Sequence:
GGG

julia> seqmatrix(seqs, :site, UInt8)
4×3 Array{UInt8,2}:
 0x01  0x01  0x01
 0x08  0x08  0x08
 0x02  0x02  0x02
 0x04  0x04  0x04

julia> seqmatrix(seqs, :seq, UInt8)
3×4 Array{UInt8,2}:
 0x01  0x08  0x02  0x04
 0x01  0x08  0x02  0x04
 0x01  0x08  0x02  0x04
```
"""
function seqmatrix(::Type{T}, vseq::AbstractVector{LongSequence{A}}, major::Symbol) where {T,A<:Alphabet}
    nseqs = length(vseq)
    @assert nseqs > 0 throw(ArgumentError("Vector of BioSequence{$A} is empty."))
    nsites = length(vseq[1])
    @inbounds for i in 2:nseqs
        length(vseq[i]) == nsites || throw(ArgumentError("Sequences in vseq must be of same length."))
    end
    if major == :site
        mat = Matrix{T}(undef, (nseqs, nsites))
        @inbounds for seq in 1:nseqs, site in 1:nsites
            mat[seq, site] = convert(T, reinterpret(UInt8, vseq[seq][site]))
        end
        return mat
    elseif major == :seq
        mat = Matrix{T}(undef, (nsites, nseqs))
        @inbounds for seq in 1:nseqs, site in 1:nsites
            mat[site, seq] = convert(T, reinterpret(UInt8, vseq[seq][site]))
        end
        return mat
    else
        throw(ArgumentError("major must be :site or :seq"))
    end
end

# Consensus
# ---------

"""
    majorityvote(seqs::AbstractVector{LongSequence{A}}) where {A<:NucleicAcidAlphabet}

Construct a sequence that is a consensus of a vector of sequences.

The consensus is established by a simple majority vote rule, where ambiguous
nucleotides cast an equal vote for each of their possible states.
For each site a winner(s) out of A, T(U), C, or G is determined, in the cases
of ties the ambiguity symbol that unifies all the winners is returned.
E.g if A and T tie, then W is inserted in the consensus. If all A, T, C, and G
tie at a site, then N is inserted in the consensus. Note this means that if a
nucletide e.g. 'C' and a gap '-' draw, the nucleotide will always win over the
gap, even though they tied.

# Examples
```julia
julia> seqs = [dna"CTCGATCGATCC", dna"CTCGAAAAATCA", dna"ATCGAAAAATCG", dna"ATCGGGGGATCG"]

4-element Array{BioSequences.BioSequence{BioSequences.DNAAlphabet{4}},1}:
 CTCGATCGATCC
 CTCGAAAAATCA
 ATCGAAAAATCG
 ATCGGGGGATCG

julia> majorityvote(seqs)
12nt DNA Sequence:
MTCGAAARATCG
```
"""
function majorityvote(seqs::AbstractVector{LongSequence{A}}) where {A<:NucleicAcidAlphabet}
    mat = seqmatrix(UInt8, seqs, :site)
    nsites = size(mat, 2)
    nseqs = size(mat, 1)
    result = BioSequence{A}(nsites)
    votes = Array{Int}(undef, 16)
    @inbounds for site in 1:nsites
        fill!(votes, 0)
        for seq in 1:nseqs
            nuc = mat[seq, site]
            votes[1] += nuc == 0x00
            votes[2] += (nuc & 0x01) != 0x00
            votes[3] += (nuc & 0x02) != 0x00
            votes[5] += (nuc & 0x04) != 0x00
            votes[9] += (nuc & 0x08) != 0x00
        end
        m = maximum(votes)
        merged = 0x00
        for i in 0x01:0x10
            merged |= ifelse(votes[i] == m, i - 0x01, 0x00)
        end
        result[site] = reinterpret(eltype(A), merged)
    end
    return result
end

### Comparisons
function Base.:(==)(seq1::SeqOrView{A}, seq2::SeqOrView{A}) where {A <: Alphabet}
    length(seq1) == length(seq2) || return false

    # If they share the same data
    isempty(seq1) && return true
    if seq1.data === seq2.data && firstbitindex(seq1) == firstbitindex(seq2)
        return true
    end

    # Check all filled UInts
    (it, (ch1, ch2, rm)) = iter_chunks(seq1, seq2)
    chunks = length(it)
    state = first_state(it)
    unroll = 8
    while chunks ≥ unroll
        same = true
        for _ in 1:unroll
            ((i, j), state) = iter_inbounds(it, state)
            same &= i == j 
        end
        same || return false
        chunks -= unroll
    end
    itval = iterate(it, state)
    while itval !== nothing
        ((i, j), state) = itval
        i == j || return false
        itval = iterate(it, state)
    end

    # Check last coding UInt (or compare two zeros, if none)
    mask = UInt64(1) << (rm & 63) - 1
    return (ch1 & mask) == (ch2 & mask)
end

function Base.:(==)(seq1::LongSequence{A}, seq2::LongSequence{A}) where {A <: Alphabet}
    length(seq1) == length(seq2) || return false
    isempty(seq1) && return true
    (data1, data2) = (seq1.data, seq2.data)

    # Check all filled UInts
    nextind = nextposition(lastbitindex(seq1))
    last_chunk_index = index(nextind) % Int - 1
    i = 1
    unroll = 8
    while i + unroll - 2 < last_chunk_index
        same = true
        @inbounds for j in 0:unroll-1
            same &= data1[i + j] == data2[i + j]
        end
        same || return false
        i += unroll
    end
    @inbounds for i in i:last_chunk_index
        data1[i] == data2[i] || return false
    end
    # Check last coding UInt, if any
    @inbounds if !iszero(offset(nextind))
        mask = UInt64(1) << (offset(nextind) & 63) - 1
        i = index(nextind)
        (data1[i] & mask) == (data2[i] & mask) || return false
    end
    return true
end

## Search
# We only dispatch on known alphabets, because new alphabets may implement == in surprising ways
function _findnext(
    cmp::Base.Fix2{<:Union{typeof(==), typeof(isequal)}},
    seq::SeqOrView{<:KNOWN_ALPHABETS},
    i::Int,
)
    enc = tryencode(Alphabet(seq), cmp.x)
    enc === nothing && return nothing
    u_enc = enc * encoding_expansion(BitsPerSymbol(seq))
    f = x -> set_zero_encoding(BitsPerSymbol(seq), x ⊻ u_enc)
    _findnext_nonzero(f, seq, i)
end

function _findnext(
    cmp::Union{<:Base.Fix2{typeof(!=)}, ComposedFunction{typeof(!), <:Base.Fix2{typeof(isequal)}}},
    seq::SeqOrView{<:KNOWN_ALPHABETS},
    i::Int,
)
    enc = tryencode(Alphabet(seq), cmp.x)
    enc === nothing && return nothing
    u_enc = enc * encoding_expansion(BitsPerSymbol(seq))
    f = x -> x ⊻ u_enc
    _findnext_nonzero(f, seq, i)
end

@inline function _findnext_nonzero(
    f,
    seq::SeqOrView{<:KNOWN_ALPHABETS},
    i::Int,
)
    vw = @inbounds view(seq, i:lastindex(seq))
    res = _findfirst_nonzero(f, vw)
    res === nothing ? nothing : res + i - 1
end


# Find first index of seq, where f is a function that takes an encodung UInt64
# and will zero out every coding element we are not looking for.
function _findfirst_nonzero(f::Function, seq::SeqOrView{<:KNOWN_ALPHABETS})
    data = seq.data
    ((head, head_bits), (body_i, body_stop), (tail, tail_bits)) = parts(seq)
    bps = bits_per_symbol(Alphabet(seq))
    
    # Use f to zero out all coding bits in head, and then mask to make sure all
    # noncoding bits are zero. Then, the chunk will be all zero if the element
    # we are not looking for is here.
    mask = (UInt(1) << (head_bits & 63)) - 1
    tz = trailing_zeros(f(head) & mask)
    tz < 64 && return div(tz % UInt, bps) % Int + 1
    # i is the current position. We just checked the head.
    i = div(head_bits % UInt, bps) % Int + 1

    # Now, process 8*64 bits at a time to improve SIMD, simply
    # skipping them if they are all zero.
    unroll = 8
    while body_i + unroll - 2 < body_stop
        u = zero(UInt64)
        for j in 0:unroll - 1
            u |= f(@inbounds(data[body_i + j]))
        end
        iszero(u) || break
        i += symbols_per_data_element(seq) * unroll
        body_i += unroll
    end

    # Check the remaining full chunks
    while body_i ≤ body_stop
        u = f(@inbounds(data[body_i]))
        if !iszero(u)
            tz = trailing_zeros(u)
            return i + div(tz % UInt, bps) % Int
        end
        body_i += 1
        i += symbols_per_data_element(seq)
    end
    
    # This is similar to the head
    mask = (UInt(1) << (tail_bits & 63)) - 1
    tz = trailing_zeros(f(tail) & mask)
    tz < 64 && return div(tz % UInt, bps) % Int + i
    nothing
end

function _findprev(
    cmp::Base.Fix2{<:Union{typeof(==), typeof(isequal)}},
    seq::SeqOrView{<:KNOWN_ALPHABETS},
    i::Int,
)
    enc = tryencode(Alphabet(seq), cmp.x)
    enc === nothing && return nothing
    u_enc = enc * encoding_expansion(BitsPerSymbol(seq))
    f = x -> set_zero_encoding(BitsPerSymbol(seq), x ⊻ u_enc)
    _findprev_nonzero(f, seq, i)
end

function _findprev(
    cmp::Union{<:Base.Fix2{typeof(!=)}, ComposedFunction{typeof(!), <:Base.Fix2{typeof(isequal)}}},
    seq::SeqOrView{<:KNOWN_ALPHABETS},
    i::Int,
)
    enc = tryencode(Alphabet(seq), cmp.x)
    enc === nothing && return nothing
    u_enc = enc * encoding_expansion(BitsPerSymbol(seq))
    f = x -> x ⊻ u_enc
    _findprev_nonzero(f, seq, i)
end

function _findprev_nonzero(
    f,
    seq::SeqOrView{<:KNOWN_ALPHABETS},
    i::Int,
)
    vw = @inbounds view(seq, 1:i)
    _findlast_nonzero(f, vw)
end

# See comments in the findfirst version of this function
function _findlast_nonzero(f::Function, seq::SeqOrView{<:KNOWN_ALPHABETS})
    data = seq.data
    ((head, head_bits), (body_stop, body_i), (tail, tail_bits)) = parts(seq)
    i = lastindex(seq)
    bps = bits_per_symbol(Alphabet(seq))
    # This part is slightly different, because the coding bits are shifted to the right,
    # but we need to count the leading bits.
    # So, we set all the unused bits to zero, then count leading zeros, and then
    # ignore the unused zero bits.
    mask = (UInt(1) << (tail_bits & 63)) - 1
    lz = leading_zeros(f(tail) & mask)
    if lz < 64
        # Compensate for@inbounds(data[body_i]) noncoding zero bits
        zero_symbols = div((lz - (64 - tail_bits)) % UInt, bps) % Int
        return i - zero_symbols
    end
    i -= div(tail_bits % UInt, bps) % Int
    unroll = 8
    # Cast to signed to prevent overflow
    (body_i, body_stop) = (body_i % Int, body_stop % Int)
    while body_i - unroll + 2 > body_stop
        u = zero(UInt64)
        for j in 0:unroll - 1
            u |= f(@inbounds(data[body_i - j]))
        end
        iszero(u) || break
        i -= symbols_per_data_element(seq) * unroll
        body_i -= unroll
    end
    while body_i ≥ body_stop
        u = f(@inbounds(data[body_i]))
        if !iszero(u)
           lz = leading_zeros(u)
           leading_symbols = (div(lz % UInt, bps) % Int)
           return i - leading_symbols
        end
        body_i -= 1
        i -= symbols_per_data_element(seq)
    end
    mask = (UInt(1) << (head_bits & 63)) - 1
    lz = leading_zeros(f(head) & mask)
    if lz < 64
        zero_symbols = div((lz - (64 - head_bits)) % UInt, bps) % Int
        return i - zero_symbols
    end
    nothing
end

encoding_expansion(::BitsPerSymbol{8}) = 0x0101010101010101
encoding_expansion(::BitsPerSymbol{4}) = 0x1111111111111111
encoding_expansion(::BitsPerSymbol{2}) = 0x5555555555555555

# For every 8-bit chunk, if the chunk is all zeros, set the lowest bit in the chunk,
# else, zero the chunk.
# E.g. 0x_0a_b0_0c_00_fe_00_ff_4e -> 0x_00_00_00_01_00_01_00_00
function set_zero_encoding(B::BitsPerSymbol{8}, enc::UInt64)
    enc = ~enc
    enc &= enc >> 4
    enc &= enc >> 2
    enc &= enc >> 1
    enc & encoding_expansion(B)
end

function set_zero_encoding(B::BitsPerSymbol{4}, enc::UInt64)
    enc = ~enc
    enc &= enc >> 2
    enc &= enc >> 1
    enc & encoding_expansion(B)
end

function set_zero_encoding(B::BitsPerSymbol{2}, enc::UInt64)
    enc = ~enc
    enc &= enc >> 1
    enc & encoding_expansion(B)
end
