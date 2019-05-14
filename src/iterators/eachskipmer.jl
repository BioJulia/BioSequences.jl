# Skipmer Iterators
# =================
#
# Iterator over all Skipmers in a biological sequence.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

eachcanonical(::Type{Skipmer{U,A,M,N,K}}, seq::BioSequence) where {U,A,M,N,K} = CanonicalSkipmers(Skipmer{U,A,M,N,K}, seq) 

struct CanonicalSkipmers{SK<:Skipmer,U<:Unsigned,S<:BioSequence}
    seq::S
    cycle_pos::Vector{UInt8}
    last_unknown::Vector{Int64}
    fkmer::Vector{U}
    rkmer::Vector{U}
    
    function CanonicalSkipmers(::Type{SK}, seq::NucleotideSeq) where {SK<:Skipmer}
        checkskipmer(SK)
        if eltype(seq) ∉ (DNA, RNA)
            throw(ArgumentError("element type must be either DNA or RNA nucleotide"))
        end
        if span(SK) > length(seq)
            throw(ArgumentError(string("The span of ", SK, " (", span(SK), ") is greater than the input sequence length (", length(seq), ").")))
        end
        N = cycle_len(SK)
        # Storage for the next N skipers being built simultaneously by the iterator.
        fkmer = Vector{encoded_data_eltype(SK)}(undef, N)
        rkmer = Vector{encoded_data_eltype(SK)}(undef, N)
        # Vector to keep track of position of last ambiguous nucleotide we saw,
        # for each of the N skipmers being built simultaneously by the iterator.
        last_unknown = Vector{Int64}(undef, N)
        # For each of the next N skipmers being build simultaneously by the iterator,
        # keep track of a cycle_pos variable: it keeps track of which nucleotides
        # in the sequence are appended to said skipmer, and is used by the
        # _consider_position! method.
        cycle_pos = Vector{UInt8}(undef, N)
        return new{SK,encoded_data_eltype(SK),typeof(seq)}(seq, cycle_pos, last_unknown, fkmer, rkmer)
    end
end

@inline Base.IteratorSize(::Type{T}) where T <: CanonicalSkipmers = Base.HasLength()
@inline Base.IteratorEltype(::Type{T}) where T <: CanonicalSkipmers = Base.HasEltype()
@inline Base.eltype(::Type{CanonicalSkipmers{SK,U,SQ}}) where {SK<:Skipmer,U<:Unsigned,SQ<:BioSequence} = SK
@inline Base.length(it::CanonicalSkipmers) = length(it.seq) - span(eltype(it)) + 1

@inline kmersize(x::CanonicalSkipmers) = kmersize(eltype(x))

@inline function kmer_mask(x::CanonicalSkipmers{SK,UT}) where {SK <: Skipmer, UT <: Unsigned}
    return (UT(1) << (kmersize(SK) * 2)) - 1
end

@inline function init_iterator!(it::CanonicalSkipmers)
    N = cycle_len(eltype(it))
    @inbounds for i in 1:N
        it.cycle_pos[i] = N - i
        it.last_unknown[i] = -1
        it.fkmer[i] = 0
        it.rkmer[i] = 0
    end
end

# Iterating over 2 bit nucleotide sequences
# -----------------------------------------

function Base.iterate(it::CanonicalSkipmers{SK, UT, SQ}) where
        {SK, UT, SQ <: BioSequence{<:NucleicAcidAlphabet{2}}}
    # This method is run the first time the iterator is called in the loop.
    
    # The length (in bp) of the cycles used to construct the skipmers.
    # Because skipmers overlap, the iterator builds the next N skipers
    # simultaneously as it passes over each base of the reference sequence.
    N = cycle_len(SK)
    # The number of bases per cycle included in very skipmer.
    M = bases_per_cycle(SK)
    # The span of the skipmer 
    S = span(SK)
    
    init_iterator!(it)
    # Take the first S bases of the sequence and "consider them".
    # The _consider_position! method takes a base at a position, it decides
    # which of the N skipmers being built should get that base, and then it
    # appends the base to those skipmers. 
    for pos in 1:S
        _consider_position!(it, pos)
    end
    # After passing over the first S bases in the sequence. Now the first of the
    # N skipers the iterator is building is now done. So now we canonicalise the
    # skipmer (whichever is less: the skipmer or its reverse complement) and
    # return it.
    fkmer = first(it.fkmer)
    rkmer = first(it.rkmer)
    outkmer = ifelse(fkmer < rkmer, fkmer, rkmer)
    return SK(outkmer), (S + 1, UInt(1))
end

function Base.iterate(it::CanonicalSkipmers{SK, UT, SQ}, state::Tuple{UInt, UInt}) where
        {SK, UT, SQ <: BioSequence{<:NucleicAcidAlphabet{2}}}
    # This method is run after the first time the iterator is called in the loop.
        
    # The length (in bp) of the cycles used to construct the skipmers.
    # Because skipmers overlap, the iterator builds the next N skipers
    # simultaneously as it passes over each base of the reference sequence.
    N = cycle_len(SK)
    # The number of bases per cycle included in very skipmer.
    M = bases_per_cycle(SK)
    
    pos = state[1]
    fi  = state[2]
    
    # If you've hit the end of the sequence, we're done.
    if pos > lastindex(it.seq)
        return nothing
    end
    # _consider_position! takes the base at pos, it decides
    # which of the N skipmers being built should get that base, and then it
    # appends the base to those skipmers.
    _consider_position!(it, pos)
    # Now we take the next skipmer that will now be complete. Get its canonical
    # form, and then return it.
    fi += 1
    fi = ifelse(fi == (N + 1), UInt(1), fi)
    fkmer = it.fkmer[fi]
    rkmer = it.rkmer[fi]
    outkmer = ifelse(fkmer < rkmer, fkmer, rkmer)
    return SK(outkmer), (pos + 1, fi)  
end

@inline function _consider_position!(it::CanonicalSkipmers{SK, UT, SQ}, pos) where
        {SK, UT, SQ <: BioSequence{<:NucleicAcidAlphabet{2}}}
    
    # The length (in bp) of the cycles used to construct the skipmers.
    # Because skipmers overlap, the iterator builds the next N skipers
    # simultaneously as it passes over each base of the reference sequence.
    N = cycle_len(SK)
    # The number of bases per cycle included in very skipmer.
    M = bases_per_cycle(SK)
    
    # For each of the N next skipmers we are building...
    for ni in 1:N
        # Update the ni'th skipmer's cycle_pos value. This is used to decide if
        # the ni'th skipmer will get the base at the position currently being considered.
        it.cycle_pos[ni] += 1
        if it.cycle_pos[ni] == N
            it.cycle_pos[ni] = 0
        end
        # If the ni'th skipmer should get the base at the position currently being
        # considered, then append said base to said skipmer. We do this using some
        # bit-flipping. We also build both the forward and reverse forms of the N
        # skipmers, at the same time for efficiency.
        if it.cycle_pos[ni] < M
            fbits = extract_encoded_element(bitindex(it.seq, pos), encoded_data(it.seq))
            _append_skipmer_bits!(it, fbits, ni)
        end
    end
end

@inline function _append_skipmer_bits!(it::CanonicalSkipmers{SK, UT, SQ}, fbits::Unsigned, cycle) where {SK, UT, SQ}
    rbits = ~fbits & typeof(fbits)(0x03)
    @inbounds it.fkmer[cycle] = ((it.fkmer[cycle] << 0x02) | fbits) & kmer_mask(it)
    @inbounds it.rkmer[cycle] = (it.rkmer[cycle] >> 0x02) | (UT(rbits) << unsigned(offset(SK, 1)))
end

# Iterating over 4 bit nucleotide sequences
# -----------------------------------------

function Base.iterate(it::CanonicalSkipmers{SK, UT, SQ}) where 
        {SK, UT, SQ <: BioSequence{<:NucleicAcidAlphabet{4}}}
    init_iterator!(it)
    pos = firstindex(it.seq)
    fi = 0x00
    return _iterate_kernel!(it, pos, fi)
end

function Base.iterate(it::CanonicalSkipmers{SK, UT, SQ}, state) where
        {SK, UT, SQ <: BioSequence{<:NucleicAcidAlphabet{4}}}
    return _iterate_kernel!(it, state[1], state[2])
end

function _iterate_kernel!(it::CanonicalSkipmers{SK, UT, SQ}, pos, fi) where
        {SK, UT, SQ <: BioSequence{<:NucleicAcidAlphabet{4}}}
    S = span(SK)
    N = cycle_len(SK)
    lastpos = lastindex(it.seq)
    # For every position in the sequence...
    while pos <= lastpos
        # Consider the base at this position, and append it to the relevant 
        # skipmers of the N next skipers being built by the iterator.
        _consider_position!(it, pos)
            
        # If we are at pos, the skip-mer that started at pos-S is now done... 
        if pos >= S
            fi += 0x01
            if fi == (N + 1)
                fi = 0x01
            end
            # If the currently finished skipmer, does not contain an ambiguous
            # nucleotide, then get the canonical form of the skipmer and return it.
            if it.last_unknown[fi] + S <= pos
                fkmer = it.fkmer[fi]
                rkmer = it.rkmer[fi]
                outkmer = ifelse(fkmer <= rkmer, fkmer, rkmer)
                return SK(outkmer), (pos + 1, fi)
            end
        end
            
        pos += 1
            
    end
    # pos > lastpos, so we are done iterating through the sequence.    
    return nothing
end

@inline function _consider_position!(it::CanonicalSkipmers{SK, UT, SQ}, pos) where
        {SK, UT, SQ <: BioSequence{<:NucleicAcidAlphabet{4}}}
        
    # The length (in bp) of the cycles used to construct the skipmers.
    # Because skipmers overlap, the iterator builds the next N skipers
    # simultaneously as it passes over each base of the reference sequence.
    N = cycle_len(SK)
    # The number of bases per cycle included in very skipmer.
    M = bases_per_cycle(SK)
    
    # For each of the N next skipmers we are building...
    for ni in 1:N
        # Update the ni'th skipmer's cycle_pos value. This is used to decide if
        # the ni'th skipmer will get the base at the position currently being considered.
        it.cycle_pos[ni] += 1
        if it.cycle_pos[ni] == N
            it.cycle_pos[ni] = 0
        end
            
        # If the ni'th skipmer should get the base at the position currently being
        # considered, then append said base to said skipmer. We do this using some
        # bit-flipping. We also build both the forward and reverse forms of the N
        # skipmers, at the same time for efficiency.
        
        # We also check that the current site is not an ambiguous nucleotide, which
        # can sometimes happen for any T<:BioSequence{<:NucleicAcidAlphabet{4}}.
        # If the base at the current position is ambigous, we update "last_unknown"
        # to mark it. Skipmers covering this position will be passed over by the
        # iterator. 
        if it.cycle_pos[ni] < M
            fbits = BioSequences.twobitnucs[extract_encoded_element(bitindex(it.seq, pos), encoded_data(it.seq)) + 0x01]
            if fbits == 0xFF
                it.last_unknown[ni] = pos
                fbits = 0x00
            end
            _append_skipmer_bits!(it, fbits, ni)
        end
    end
end
