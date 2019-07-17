### Skipmer Iterators
###
###
### Iterator over all Skipmers in a biological sequence.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

struct SkipmerFactoryResult{T<:AbstractMer}
    position::Int
    fw::T
    bw::T
end

mutable struct SkipmerFactory{S<:LongNucleotideSequence,U<:Unsigned,K}
    seq::S
    cycle_pos::Vector{UInt8}
    last_unknown::Vector{Int64}
    fkmer::Vector{U}
    rkmer::Vector{U}
    position::Int
    finished::Int
    n::Int
    cycle_len::Int
    bases_per_cycle::Int
    span::Int
    
    function SkipmerFactory(::Type{T},
                            seq::LongSequence{A},
                            bases_per_cycle::Int = 2,
                            cycle_len::Int = 3) where {A<:NucleicAcidAlphabet,T<:AbstractMer}
        
        checkmer(T)
        
        span = UInt(ceil(cycle_len * (ksize(T) / bases_per_cycle - 1) + bases_per_cycle))
        
        if eltype(seq) ∉ (DNA, RNA)
            throw(ArgumentError("element type must be either DNA or RNA nucleotide"))
        elseif bases_per_cycle > cycle_len
            throw(ArgumentError("bases per cycle must not be greater than the cycle length"))
        elseif span > length(seq)
            throw(ArgumentError(string("span of each skipmer (", span, ") is greater than the input sequence length (", length(seq), ").")))
        elseif eltype(seq) != eltype(T)
            throw(ArgumentError(string("skipmer type chosen must have same element type as chosen sequence")))
        end
        
        # For each of the next N skipmers being build simultaneously by the iterator,
        # keep track of a cycle_pos variable: it keeps track of which nucleotides
        # in the sequence are appended to said skipmer, and is used by the
        # _consider_position! method.
        cycle_pos = Vector{UInt8}(undef, cycle_len)
        
        # Vector to keep track of position of last ambiguous nucleotide we saw,
        # for each of the N skipmers being built simultaneously by the iterator.
        last_unknown = Vector{Int64}(undef, cycle_len)
        
        # Storage for the next N skipers being built simultaneously by the iterator.
        U = encoded_data_type(T)
        fkmer = Vector{U}(undef, cycle_len)
        rkmer = Vector{U}(undef, cycle_len)
        
        gen = new{LongSequence{A},U,ksize(T)}(seq, cycle_pos, last_unknown, fkmer, rkmer, 0, 0, 0, cycle_len, bases_per_cycle, span)
        _init_generator!(gen)
        
        return gen
    end
end

@inline function mertype(::Type{SkipmerFactory{S,UInt64,K}}) where {S<:LongNucleotideSequence,K}
    return Mer{minimal_alphabet(typeof(Alphabet(S))),K}
end

@inline function mertype(::Type{SkipmerFactory{S,UInt128,K}}) where {S<:LongNucleotideSequence,K}
    return BigMer{minimal_alphabet(typeof(Alphabet(S))),K}
end

@inline mertype(gen::SkipmerFactory) = mertype(typeof(gen))

@inline function Base.eltype(::Type{S}) where {S<:SkipmerFactory}
    return SkipmerFactoryResult{mertype(S)}
end



@inline Base.IteratorSize(::Type{T}) where T <: SkipmerFactory = Base.HasLength()

@inline Base.IteratorEltype(::Type{T}) where T <: SkipmerFactory = Base.HasEltype()

@inline Base.length(gen::SkipmerFactory)::Int = length(gen.seq) - gen.span + 1

function set_sequence!(gen::SkipmerFactory{S}, seq::S) where {S<:LongNucleotideSequence}
    gen.seq = seq
    _init_generator!(gen)
end

function _init_generator!(gen::SkipmerFactory)
    N = gen.cycle_len
    
    @inbounds for i in 1:N
        gen.cycle_pos[i] = N - i
        gen.last_unknown[i] = -1
        gen.fkmer[i] = 0
        gen.rkmer[i] = 0
    end
    
    gen.finished = 0
    gen.position = 0
    
    for i in 1:gen.span-1
        _consider_next_position!(gen)
    end
end

@inline function _get_base_bits(seq::LongSequence{<:NucleicAcidAlphabet{2}}, position)
    return extract_encoded_element(bitindex(seq, position), encoded_data(seq))
end

@inline function _get_base_bits(seq::LongSequence{<:NucleicAcidAlphabet{4}}, position)
    return twobitnucs[extract_encoded_element(bitindex(seq, position), encoded_data(seq)) + 0x01]
end

@inline function _append_mers!(gen::SkipmerFactory{LongSequence{A},U,K}, fbits::U, rbits::U) where {A<:NucleicAcidAlphabet{2},U,K}
    # For each of the N next skipmers we are building...
    N = gen.cycle_len
    M = gen.bases_per_cycle
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

@inline function _append_mers!(gen::SkipmerFactory{LongSequence{A},U,K}, fbits::U, rbits::U) where {A<:NucleicAcidAlphabet{4},U,K}
    N = gen.cycle_len
    M = gen.bases_per_cycle
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

@inline function _consider_next_position!(gen::SkipmerFactory{S,U,K}) where {S,U,K}
    nextposition = gen.position + 1
    fbits = U(_get_base_bits(gen.seq, nextposition))
    gen.position = nextposition
    rbits = ~fbits & U(0x03)
    _append_mers!(gen, fbits, rbits)
end

function nextmer(gen::SkipmerFactory{LongSequence{A},U,K}) where {A<:NucleicAcidAlphabet{2},U,K}
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
    nextfinished = ifelse(nextfinished > gen.cycle_len, 1, nextfinished)
    gen.finished = nextfinished
    @inbounds fkmer = gen.fkmer[nextfinished]
    @inbounds rkmer = gen.rkmer[nextfinished]
    newn = gen.n + 1
    gen.n = newn
    return SkipmerFactoryResult(newn, mertype(gen)(fkmer), mertype(gen)(rkmer))
end

function nextmer(gen::SkipmerFactory{LongSequence{A},U,K}) where {A<:NucleicAcidAlphabet{4},U,K}
    lastpos = lastindex(gen.seq)
    span = gen.span
    # For every position in the sequence...
    while gen.position < lastpos
        # Consider the base at the next position, and append it to the relevant 
        # skipmers of the N next skipers being built by the iterator.
        _consider_next_position!(gen)
            
        # If we are at pos, the skip-mer that started at pos-S is now done... 
        if gen.position >= span
            nextfinished = gen.finished + 1
            gen.finished = ifelse(nextfinished > gen.cycle_len, 1, nextfinished)
            # If the currently finished skipmer, does not contain an ambiguous
            # nucleotide, then get the canonical form of the skipmer and return it.
            newn = gen.n + 1
            gen.n = newn
            if gen.last_unknown[gen.finished] + span <= gen.position
                fkmer = gen.fkmer[gen.finished]
                rkmer = gen.rkmer[gen.finished]
                return SkipmerFactoryResult(newn, mertype(gen)(fkmer), mertype(gen)(rkmer))
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

@inline function each(::Type{T}, seq::BioSequence, bases_per_cycle::Int, cycle_len::Int) where {T<:AbstractMer}
    return SkipmerFactory(T, seq, bases_per_cycle, cycle_len)
end
