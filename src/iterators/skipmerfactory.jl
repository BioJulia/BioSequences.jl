# Skipmer Iterators
# =================
#
# Iterator over all Skipmers in a biological sequence.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

mutable struct SkipmerFactory{S<:LongNucleotideSequence,U<:Unsigned,M,N,K}
    seq::S
    cycle_pos::Vector{UInt8}
    last_unknown::Vector{Int64}
    fkmer::Vector{U}
    rkmer::Vector{U}
    position::Int
    finished::Int
    n::Int
    
    function SkipmerFactory(::Type{Skipmer{U,A,M,N,K}}, seq::S) where {S<:LongNucleotideSequence,U<:Unsigned,A,M,N,K}
        checkskipmer(Skipmer{U,A,M,N,K})
        
        if _span(M,N,K) > length(seq)
            throw(ArgumentError(string("The span of each skipmer (", _span(M,N,K), ") is greater than the input sequence length (", length(seq), ").")))
        end
        
        if eltype(Skipmer{U,A,M,N,K}) !== eltype(seq)
            throw(ArgumentError(string("Input is a ", eltype(seq), "sequence, can't generate ", eltype(Skipmer{U,A,M,N,K}), "skipmers.")))
        end
        
        # For each of the next N skipmers being build simultaneously by the iterator,
        # keep track of a cycle_pos variable: it keeps track of which nucleotides
        # in the sequence are appended to said skipmer, and is used by the
        # _consider_position! method.
        cycle_pos = Vector{UInt8}(undef, N)
        
        # Vector to keep track of position of last ambiguous nucleotide we saw,
        # for each of the N skipmers being built simultaneously by the iterator.
        last_unknown = Vector{Int64}(undef, N)
        
        # Storage for the next N skipers being built simultaneously by the iterator.
        fkmer = Vector{U}(undef, N)
        rkmer = Vector{U}(undef, N)
        
        gen = new{S,U,M,N,K}(seq, cycle_pos, last_unknown, fkmer, rkmer, 0, 0, 0)
        _init_generator!(gen)
        
        return gen
    end
end

@inline function mertype(::Type{SkipmerFactory{S,U,M,N,K}}) where {U<:Unsigned,M,N,K,S<:LongNucleotideSequence}
    return Skipmer{U,ifelse(eltype(S) === DNA, DNAAlphabet{2}, RNAAlphabet{2}),M,N,K}
end
@inline mertype(gen::SkipmerFactory) = mertype(typeof(gen))

@inline function Base.eltype(::Type{SkipmerFactory{S,U,M,N,K}}) where {U<:Unsigned,M,N,K,S<:LongNucleotideSequence}
    ST = mertype(SkipmerFactory{S,U,M,N,K})
    return Tuple{Int,ST,ST}
end

@inline Base.IteratorSize(::Type{T}) where T <: SkipmerFactory = Base.HasLength()

@inline Base.IteratorEltype(::Type{T}) where T <: SkipmerFactory = Base.HasEltype()

@inline Base.length(gen::SkipmerFactory)::Int = length(gen.seq) - span(mertype(gen)) + 1

function set_sequence!(gen::SkipmerFactory{S}, seq::S) where {S<:LongNucleotideSequence}
    gen.seq = seq
    _init_generator!(gen)
end

function _init_generator!(gen::SkipmerFactory)
    N = cycle_len(mertype(gen))
    
    @inbounds for i in 1:N
        gen.cycle_pos[i] = N - i
        gen.last_unknown[i] = -1
        gen.fkmer[i] = 0
        gen.rkmer[i] = 0
    end
    
    gen.finished = 0
    gen.position = 0
    
    for i in 1:span(mertype(gen))-1
        _consider_next_position!(gen)
    end
end

@inline function _get_base_bits(seq::LongSequence{<:NucleicAcidAlphabet{2}}, position)
    return extract_encoded_element(bitindex(seq, position), encoded_data(seq))
end

@inline function _get_base_bits(seq::LongSequence{<:NucleicAcidAlphabet{4}}, position)
    return twobitnucs[extract_encoded_element(bitindex(seq, position), encoded_data(seq)) + 0x01]
end

@inline function _append_mers!(gen::SkipmerFactory{LongSequence{A},U,M,N,K}, fbits::U, rbits::U) where {A<:NucleicAcidAlphabet{2},U,M,N,K}
    # For each of the N next skipmers we are building...
    for ni in 1:N
        # Update the ni'th skipmer's cycle_pos value.
        # This is used to decide if the ni'th skipmer will get the base at the
        # position currently being considered.
        nextcycle = gen.cycle_pos[ni] + 1
        gen.cycle_pos[ni] = ifelse(nextcycle == N, 0, nextcycle)
        # If the ni'th skipmer should get the base at the position currently being
        # considered, then append said base to said skipmer. We do this using some
        # bit-flipping. We also build both the forward and reverse forms of the N
        # skipmers, at the same time for efficiency.
        if gen.cycle_pos[ni] < M
            gen.fkmer[ni] = ((gen.fkmer[ni] << 0x02) | fbits)
            gen.rkmer[ni] = (gen.rkmer[ni] >> 0x02) | (rbits << 2(K - 1))
        end
    end
end

@inline function _append_mers!(gen::SkipmerFactory{LongSequence{A},U,M,N,K}, fbits::U, rbits::U) where {A<:NucleicAcidAlphabet{4},U,M,N,K}
    if fbits == 0xFF
        for ni in 1:N
            nextcycle = gen.cycle_pos[ni] + 1
            gen.cycle_pos[ni] = ifelse(nextcycle == N, 0, nextcycle)
            if gen.cycle_pos[ni] < M
                gen.fkmer[ni] = (gen.fkmer[ni] << 0x02)
                gen.rkmer[ni] = (gen.rkmer[ni] >> 0x02)
            end
        end
    else
        for ni in 1:N
            nextcycle = gen.cycle_pos[ni] + 1
            gen.cycle_pos[ni] = ifelse(nextcycle == N, 0, nextcycle)
            if gen.cycle_pos[ni] < M
                gen.fkmer[ni] = ((gen.fkmer[ni] << 0x02) | fbits)
                gen.rkmer[ni] = (gen.rkmer[ni] >> 0x02) | (rbits << 2(K - 1))
            end
        end
    end
end

@inline function _consider_next_position!(gen::SkipmerFactory{S,U,M,N,K}) where {S,U,M,N,K}
    nextposition = gen.position + 1
    fbits = U(_get_base_bits(gen.seq, nextposition))
    gen.position = nextposition
    rbits = ~fbits & U(0x03)
    _append_mers!(gen, fbits, rbits)
end

function nextmer(gen::SkipmerFactory{LongSequence{A},U,M,N,K}) where {A<:NucleicAcidAlphabet{2},U,M,N,K}
    # If you've hit the end of the sequence, we're done.
    if gen.position == lastindex(gen.seq)
        return nothing
    end
    # _consider_next_position! takes the base at pos, it decides
    # which of the N skipmers being built should get that base, and then it
    # appends the base to those skipmers.
    _consider_next_position!(gen)
    # Now we take the next skipmer that will now be complete. Get its canonical
    # form, and then return it.
    nextfinished = gen.finished + 1
    nextfinished = ifelse(nextfinished > N, 1, nextfinished)
    gen.finished = nextfinished
    @inbounds fkmer = gen.fkmer[nextfinished]
    @inbounds rkmer = gen.rkmer[nextfinished]
    newn = gen.n + 1
    gen.n = newn
    return (newn, mertype(gen)(fkmer), mertype(gen)(rkmer))
end

function nextmer(gen::SkipmerFactory{LongSequence{A},U,M,N,K}) where {A<:NucleicAcidAlphabet{4},U,M,N,K}
    lastpos = lastindex(gen.seq)
    # For every position in the sequence...
    while gen.position < lastpos
        # Consider the base at the next position, and append it to the relevant 
        # skipmers of the N next skipers being built by the iterator.
        _consider_next_position!(gen)
            
        # If we are at pos, the skip-mer that started at pos-S is now done... 
        if gen.position >= _span(M, N, K)
            nextfinished = gen.finished + 1
            gen.finished = ifelse(nextfinished > N, 1, nextfinished)
            # If the currently finished skipmer, does not contain an ambiguous
            # nucleotide, then get the canonical form of the skipmer and return it.
            newn = gen.n + 1
            gen.n = newn
            if gen.last_unknown[gen.finished] + _span(M, N, K) <= gen.position
                fkmer = gen.fkmer[gen.finished]
                rkmer = gen.rkmer[gen.finished]
                return (newn, mertype(gen)(fkmer), mertype(gen)(rkmer))
            end
        end
    end
    # pos => lastpos, so we are done iterating through the sequence.    
    return nothing
end

@inline function Base.iterate(gen::SkipmerFactory)
    _init_generator!(gen)
    return iterate(gen, nothing)
end

@inline function Base.iterate(gen::SkipmerFactory, ::Nothing)
    x = nextmer(gen)
    if isnothing(x)
        return nothing
    else
        return x, nothing
    end
end

@inline function each(::Type{Skipmer{U,A,M,N,K}}, seq::BioSequence) where {U,A,M,N,K}
    return ((@inbounds it[1], @inbounds it[2]) for it in SkipmerFactory(Skipmer{U,A,M,N,K}, seq))
end

@inline function eachcanonical(::Type{Skipmer{U,A,M,N,K}}, seq::BioSequence) where {U,A,M,N,K}
    return ((@inbounds it[1], min(@inbounds it[2], @inbounds it[3])) for it in SkipmerFactory(Skipmer{U,A,M,N,K}, seq))
end
