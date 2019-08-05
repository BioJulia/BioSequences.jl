import BioSequences: SkipmerFactory, MerIterResult

function char_complement(::Type{DNA}, c::Char)
    c === 'A' && return 'T'
    c === 'T' && return 'A'
    c === 'C' && return 'G'
    c === 'G' && return 'C'
    error("Oh god, wrong character!")
end

function char_complement(::Type{RNA}, c::Char)
    c === 'A' && return 'U'
    c === 'U' && return 'A'
    c === 'C' && return 'G'
    c === 'G' && return 'C'
    error("Oh god, wrong character!")
end

function string_reverse_complement(::Type{N}, seq::AbstractString) where {N<:NucleicAcid}
    return String([char_complement(N, c) for c in reverse(seq)])
end

function string_iscanonical(::Type{N}, seq::AbstractString) where {N<:NucleicAcid}
    i = 1
    j = length(seq)
    @inbounds while i <= j
        f = seq[i]
        r = char_complement(N, seq[j])
        f < r && return true
        r < f && return false
        i += 1
        j -= 1
    end
    return true
end

function string_canonical(::Type{N}, seq::AbstractString) where {N<:NucleicAcid}
    if string_iscanonical(N, seq)
        return seq
    else
        return string_reverse_complement(N, seq)
    end
end

@testset "Kmers" begin
    function string_eachkmer(::Type{N}, seq::AbstractString, k, step) where {N<:NucleicAcid}
        kmers = Tuple{String,String}[]
        i = 1
        for i in 1:step:length(seq) - k + 1
            subseq = seq[i:i + k - 1]
            if !in('N', subseq)
                push!(kmers, (subseq, string_reverse_complement(N, subseq)))
            end
        end
        return kmers
    end
    
    function test_eachkmer(S, seq::AbstractString, k, step)
        alph = BioSequences.minimal_alphabet(Alphabet(S))
        mertype = Mer{alph, k}
        iter = each(mertype, S(seq), step)
        
        @test eltype(iter) == MerIterResult{mertype}
        
        xs = [(String(x), String(y)) for (i, x, y) in iter]
        ys = string_eachkmer(eltype(S), seq, k, step)
        
        @test xs == ys
    end
    
    function test_iteratorlength(S, seq::AbstractString, k, step)
         @test length(each(Mer{BioSequences.minimal_alphabet(Alphabet(S)),k}, S(seq), step)) == length(string_eachkmer(eltype(S), seq, k, step))
    end
    
    for k in [1, 3, 16, 32], step in 1:3, len in [1, 2, 3, 5, 10, 100, 1000]
        test_eachkmer(LongSequence{DNAAlphabet{4}}, random_dna(len), k, step)
        test_eachkmer(LongSequence{RNAAlphabet{4}}, random_rna(len), k, step)
        test_eachkmer(ReferenceSequence, random_dna(len), k, step)

        probs = [0.25, 0.25, 0.25, 0.25, 0.00]

        test_eachkmer(LongSequence{DNAAlphabet{2}}, random_dna(len, probs), k, step)
        test_eachkmer(LongSequence{RNAAlphabet{2}}, random_rna(len, probs), k, step)

        test_iteratorlength(LongSequence{DNAAlphabet{2}}, random_dna(len, probs), k, step)
        test_iteratorlength(LongSequence{RNAAlphabet{2}}, random_rna(len, probs), k, step)
    end
    
    @test isempty(collect(each(DNAMer{1}, dna"")))
    @test isempty(collect(each(DNAMer{1}, dna"NNNNNNNNNN")))
    @test_throws Exception each(DNAMer{-1}, dna"ACGT")
    @test_throws Exception each(DNAMer{0}, dna"ACGT")
    @test_throws Exception each(DNAMer{33}, dna"ACGT")
    @test collect(each(DNAMer{3}, dna"AC-TGAG--TGC")) == [MerIterResult(4, DNACodon(DNA_T, DNA_G, DNA_A), DNACodon(DNA_T, DNA_C, DNA_A)),
                                                          MerIterResult(5, DNACodon(DNA_G, DNA_A, DNA_G), DNACodon(DNA_C, DNA_T, DNA_C)),
                                                          MerIterResult(10, DNACodon(DNA_T, DNA_G, DNA_C), DNACodon(DNA_G, DNA_C, DNA_A))]
end

@testset "SkipmerFactory" begin
    seq2 = LongSequence{DNAAlphabet{2}}("GAGGGGGATGCCCCTCTTTGAGCCCAAGG")
    seq4 = LongSequence{DNAAlphabet{4}}("GAGGGGGATGCCCCTCTTTGAGCCCAAGG")
    
    #function test_iter_traits(it::CanonicalSkipmers{ST,U,S}) where {ST,U,S}
    function test_iter_traits(it::SkipmerFactory{S,UInt64,K}) where {S,K}
        A = ifelse(eltype(S) == DNA, DNAAlphabet{2}, RNAAlphabet{2})
        @test eltype(it) == MerIterResult{Mer{A,K}}
        @test Base.IteratorSize(it) == Base.HasLength()
        @test Base.IteratorEltype(it) == Base.HasEltype()
    end
    
    function test_iter_traits(it::SkipmerFactory{S,UInt128,K}) where {S,K}
        A = ifelse(eltype(S) == DNA, DNAAlphabet{2}, RNAAlphabet{2})
        @test eltype(it) == MerIterResult{BigMer{A,K}}
        @test Base.IteratorSize(it) == Base.HasLength()
        @test Base.IteratorEltype(it) == Base.HasEltype()
    end
    
    # Test when Skipmers are based on UInt64, and sample nucletides using
    # M = 2, N = 3, K = 14 pattern.
    
    #"GAGGGGGATGCCCCTCTTTGAGCCCAAGG"
    #"GA GG GA GC CC CT TG GC CA GG"
    #" AG GG AT CC CT TT GA CC AA G"
    #"G GG GG TG CC TC TT AG CC AG "
    
    fwans = ["GAGGGAGCCCCTTG", "AGGGATCCCTTTGA", "GGGGTGCCTCTTAG",
             "GGGAGCCCCTTGGC", "GGATCCCTTTGACC", "GGTGCCTCTTAGCC",
             "GAGCCCCTTGGCCA", "ATCCCTTTGACCAA", "TGCCTCTTAGCCAG",
             "GCCCCTTGGCCAGG"]
    
    M = Mer{DNAAlphabet{2},14}
    BIGM = BigMer{DNAAlphabet{2},14}
    
    ans = collect(MerIterResult(x...) for x in zip(eachindex(fwans), M.(fwans), M.(string_reverse_complement.(DNA, fwans))))
    bigans = collect(MerIterResult(x...) for x in zip(eachindex(fwans), BIGM.(fwans), BIGM.(string_reverse_complement.(DNA, fwans))))
    
    #ans = ["CAAGGGGCTCCCTC", "AGGGATCCCTTTGA", "CTAAGAGGCACCCC",
    #       "GCCAAGGGGCTCCC", "GGATCCCTTTGACC", "GGCTAAGAGGCACC",
    #       "GAGCCCCTTGGCCA", "ATCCCTTTGACCAA", "CTGGCTAAGAGGCA",
    #       "CCTGGCCAAGGGGC"]

    #test_iter_traits(eachcanonical(ST, seq2))
    test_iter_traits(SkipmerFactory(M, seq2, 2, 3))
    #test_iter_traits(eachcanonical(ST, seq4))
    test_iter_traits(SkipmerFactory(M, seq4, 2, 3))
    #@test collect(eachcanonical(ST, seq2)) == ST.(ans)
    @test collect(SkipmerFactory(M, seq2, 2, 3)) == ans
    #@test collect(eachcanonical(ST, seq4)) == ST.(ans)
    @test collect(SkipmerFactory(M, seq4, 2, 3)) == ans
    
    # Test the same pattern, but with UInt128 for big capacity skipmers.
    #test_iter_traits(eachcanonical(BIGST, seq2))
    test_iter_traits(SkipmerFactory(BIGM, seq2, 2, 3))
    #test_iter_traits(eachcanonical(BIGST, seq4))
    test_iter_traits(SkipmerFactory(BIGM, seq4, 2, 3))
    #@test collect(eachcanonical(BIGST, seq2)) == BIGST.(ans)
    @test collect(SkipmerFactory(BIGM, seq2, 2, 3)) == bigans
    #@test collect(eachcanonical(BIGST, seq4)) == BIGST.(ans)
    @test collect(SkipmerFactory(BIGM, seq4, 2, 3)) == bigans
    
    # Test when Skipmers are based on UInt64, and sample nucleotides using
    # M = 1, N = 3, K = 7 pattern.
    M2 = Mer{DNAAlphabet{2},7}
    BIGM2 = BigMer{DNAAlphabet{2},7}
    
    fwans2 = ["GGGGCCT", "AGACCTG", "GGTCTTA", "GGGCCTG",
             "GACCTGC", "GTCTTAC", "GGCCTGC", "ACCTGCA",
             "TCTTACA", "GCCTGCG", "CCTGCAG"]
    
    ans2 = collect(MerIterResult(x...) for x in zip(eachindex(fwans2), M2.(fwans2), M2.(string_reverse_complement.(DNA, fwans2))))
    bigans2 = collect(MerIterResult(x...) for x in zip(eachindex(fwans2), BIGM2.(fwans2), BIGM2.(string_reverse_complement.(DNA, fwans2))))
    
    # Does iterator stop you if the span of your skipmer exceeds the len of the
    # sequence?
    @test_throws ArgumentError SkipmerFactory(Mer{DNAAlphabet{2},14}, seq2, 1, 3)
    
    test_iter_traits(SkipmerFactory(M2, seq2, 1, 3))
    test_iter_traits(SkipmerFactory(M2, seq4, 1, 3))
    @test collect(SkipmerFactory(M2, seq2, 1, 3)) == ans2
    @test collect(SkipmerFactory(M2, seq4, 1, 3)) == ans2
    
    # Test the same pattern, but with UInt128 for big capacity skipmers.
    test_iter_traits(SkipmerFactory(BIGM2, seq2, 1, 3))
    test_iter_traits(SkipmerFactory(BIGM2, seq4, 1, 3))
    @test collect(SkipmerFactory(BIGM2, seq2, 1, 3)) == bigans2
    @test collect(SkipmerFactory(BIGM2, seq4, 1, 3)) == bigans2
    
    # Test when Skipmers are generated using M = 2, N = 4, K = 8 pattern.
    M3 = Mer{DNAAlphabet{2},8}
    BIGM3 = BigMer{DNAAlphabet{2},8}
    
    fwans3 = ["GAGGTGCC", "AGGGGCCT", "GGGACCTC", "GGATCCCT",
              "GGTGCCTT", "GGGCCTTT", "GACCTCTG", "ATCCCTGA",
              "TGCCTTAG", "GCCTTTGC", "CCTCTGCC","CCCTGACC",
              "CCTTAGCA", "CTTTGCAA", "TCTGCCAG", "CTGACCGG"]
    
    ans3 = collect(MerIterResult(x...) for x in zip(eachindex(fwans3), M3.(fwans3), M3.(string_reverse_complement.(DNA, fwans3))))
    bigans3 = collect(MerIterResult(x...) for x in zip(eachindex(fwans3), BIGM3.(fwans3), BIGM3.(string_reverse_complement.(DNA, fwans3))))
    
    test_iter_traits(SkipmerFactory(M3, seq2, 2, 4))
    test_iter_traits(SkipmerFactory(M3, seq4, 2, 4))
    @test collect(SkipmerFactory(M3, seq2, 2, 4)) == ans3
    @test collect(SkipmerFactory(M3, seq4, 2, 4)) == ans3
    
    # Test the same pattern, but with UInt128 for big capacity skipmers.
    test_iter_traits(SkipmerFactory(BIGM3, seq2, 2, 4))
    test_iter_traits(SkipmerFactory(BIGM3, seq4, 2, 4))
    @test collect(SkipmerFactory(BIGM3, seq2, 2, 4)) == bigans3
    @test collect(SkipmerFactory(BIGM3, seq4, 2, 4)) == bigans3
end

