# Basic Operators
# ---------------

function count_gc(seq::MutableBioSequence{<:Union{DNAAlphabet{2},RNAAlphabet{2}}})
    function count(x)
        # bit parallel counter of GC
        c =  x & 0x5555555555555555
        g = (x & 0xAAAAAAAAAAAAAAAA) >> 1
        return count_ones(c ⊻ g)
    end
    n = 0
    i = BitIndex(seq, 1)
    stop = BitIndex(seq, lastindex(seq) + 1)
    if offset(i) != 0 && i < stop
        # align the bit index to the beginning of a block boundary
        o = offset(i)
        n += count((seq.data[index(i)] >> o) & bitmask(stop - i))
        i += 64 - o
        @assert offset(i) == 0
    end
    while i ≤ stop - 64
        @inbounds n += count(seq.data[index(i)])
        i += 64
    end
    if i < stop
        n += count(seq.data[index(i)] & bitmask(offset(stop)))
    end
    return n
end

function count_gc(seq::MutableBioSequence{<:Union{DNAAlphabet{4},RNAAlphabet{4}}})
    function count(x)
        # bit parallel counter of GC
        a =  x & 0x1111111111111111
        c = (x & 0x2222222222222222) >> 1
        g = (x & 0x4444444444444444) >> 2
        t = (x & 0x8888888888888888) >> 3
        return count_ones((c | g) & ~(a | t))
    end
    n = 0
    i = BitIndex(seq, 1)
    stop = BitIndex(seq, lastindex(seq) + 1)
    if offset(i) != 0 && i < stop
        # align the bit index to the beginning of a block boundary
        o = offset(i)
        n += count((seq.data[index(i)] >> o) & bitmask(stop - i))
        i += 64 - o
        @assert offset(i) == 0
    end
    while i ≤ stop - 64
        @inbounds n += count(seq.data[index(i)])
        i += 64
    end
    if i < stop
        n += count(seq.data[index(i)] & bitmask(offset(stop)))
    end
    return n
end


# Site counting
# -------------

include("site_counting/site_counting.jl")


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
function seqmatrix(vseq::AbstractVector{MutableBioSequence{A}}, major::Symbol) where {A<:Alphabet}
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
function seqmatrix(::Type{T}, vseq::AbstractVector{MutableBioSequence{A}}, major::Symbol) where {T,A<:Alphabet}
    nseqs = length(vseq)
    @assert nseqs > 0 throw(ArgumentError("Vector of BioSequence{$A} is empty."))
    nsites = length(vseq[1])
    @inbounds for i in 2:nseqs
        length(vseq[i]) == nsites || throw(ArgumentError("Sequences in vseq must be of same length."))
    end
    if major == :site
        mat = Matrix{T}(undef, (nseqs, nsites))
        @inbounds for seq in 1:nseqs, site in 1:nsites
            mat[seq, site] = convert(T, vseq[seq][site])
        end
        return mat
    elseif major == :seq
        mat = Matrix{T}(undef, (nsites, nseqs))
        @inbounds for seq in 1:nseqs, site in 1:nsites
            mat[site, seq] = convert(T, vseq[seq][site])
        end
        return mat
    else
        throw(ArgumentError("major must be :site or :seq"))
    end
end

# Consensus
# ---------

"""
    majorityvote(seqs::AbstractVector{MutableBioSequence{A}}) where {A<:NucleicAcidAlphabet}

Construct a sequence that is a consensus of a vector of sequences.

The consensus is established by a simple majority vote rule, where amiguous
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
function majorityvote(seqs::AbstractVector{MutableBioSequence{A}}) where {A<:NucleicAcidAlphabet}
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
