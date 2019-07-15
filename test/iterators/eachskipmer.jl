import BioSequences: SkipmerFactory, SkipmerFactoryResult

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

function string_reverse_complement(::Type{N}, seq::AbstractString) where {N <: NucleicAcid}
    return String([char_complement(N, c) for c in reverse(seq)])
end

function string_iscanonical(::Type{N}, seq::AbstractString) where {N <: NucleicAcid}
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

function string_canonical(::Type{N}, seq::AbstractString) where {N <: NucleicAcid}
    if string_iscanonical(N, seq)
        return seq
    else
        return string_reverse_complement(N, seq)
    end
end

@testset "Kmers" begin
    function string_eachkmer(::Type{N}, seq::AbstractString, k, step, canonical = false) where {N <: NucleicAcid}
        kmers = String[]
        i = 1
        for i in 1:step:length(seq) - k + 1
            subseq = seq[i:i + k - 1]
            if !in('N', subseq)
                if canonical
                    push!(kmers, string_canonical(N, subseq))
                else
                    push!(kmers, subseq)
                end
            end
        end
        return kmers
    end
    
    function test_eachkmer(S, seq::AbstractString, k, step)
        alph = BioSequences.minimal_alphabet(Alphabet(S))
        mertype = Kmer{UInt64, alph, k}
        iter = each(mertype, S(seq), step)
        
        @test eltype(iter) == Tuple{Int,mertype}
        
        xs = [convert(String, x) for (i, x) in iter]
        ys = string_eachkmer(eltype(S), seq, k, step)
        
        @test xs == ys
    end
    
    function test_eachkmer_canonical(S, seq::AbstractString, k, step)
        alph = BioSequences.minimal_alphabet(Alphabet(S))
        mertype = Kmer{UInt64, alph, k}
        
        xs = [convert(String, x) for (i, x) in eachcanonical(mertype, S(seq), step)]
        ys = string_eachkmer(eltype(S), seq, k, step, true)
        
        @test xs == ys
    end
    
    function test_iteratorlength(S, seq::AbstractString, k, step)
         @test length(each(Kmer{UInt64,BioSequences.minimal_alphabet(Alphabet(S)),k}, S(seq), step)) == length(string_eachkmer(eltype(S), seq, k, step))
    end
    
    for k in [1, 3, 16, 32], step in 1:3, len in [1, 2, 3, 5, 10, 100, 1000]
        test_eachkmer(LongSequence{DNAAlphabet{4}}, random_dna(len), k, step)
        test_eachkmer(LongSequence{RNAAlphabet{4}}, random_rna(len), k, step)
        test_eachkmer(ReferenceSequence, random_dna(len), k, step)
        
        test_eachkmer_canonical(LongSequence{DNAAlphabet{4}}, random_dna(len), k, step)
        test_eachkmer_canonical(LongSequence{RNAAlphabet{4}}, random_rna(len), k, step)
        test_eachkmer_canonical(ReferenceSequence, random_dna(len), k, step)

        probs = [0.25, 0.25, 0.25, 0.25, 0.00]

        test_eachkmer(LongSequence{DNAAlphabet{2}}, random_dna(len, probs), k, step)
        test_eachkmer(LongSequence{RNAAlphabet{2}}, random_rna(len, probs), k, step)
        
        test_eachkmer_canonical(LongSequence{DNAAlphabet{2}}, random_dna(len, probs), k, step)
        test_eachkmer_canonical(LongSequence{RNAAlphabet{2}}, random_rna(len, probs), k, step)

        test_iteratorlength(LongSequence{DNAAlphabet{2}}, random_dna(len, probs), k, step)
        test_iteratorlength(LongSequence{RNAAlphabet{2}}, random_rna(len, probs), k, step)
    end
    
    @test isempty(collect(each(DNAKmer{1}, dna"")))
    @test isempty(collect(each(DNAKmer{1}, dna"NNNNNNNNNN")))
    @test_throws Exception each(DNAKmer{-1}, dna"ACGT")
    @test_throws Exception each(DNAKmer{0}, dna"ACGT")
    @test_throws Exception each(DNAKmer{33}, dna"ACGT")
    @test_throws Exception eachcanonical(DNAKmer{-1}, dna"ACGT")
    @test_throws Exception eachcanonical(DNAKmer{0}, dna"ACGT")
    @test_throws Exception eachcanonical(DNAKmer{33}, dna"ACGT")
    @test collect(each(DNAKmer{3}, dna"AC-TGAG--TGC")) == [(4, DNACodon(DNA_T, DNA_G, DNA_A)),
                                                           (5, DNACodon(DNA_G, DNA_A, DNA_G)),
                                                           (10, DNACodon(DNA_T, DNA_G, DNA_C))]
    
end

#@testset "CanonicalSkipmers" begin
@testset "SkipmerFactory" begin
    seq2 = LongSequence{DNAAlphabet{2}}("GAGGGGGATGCCCCTCTTTGAGCCCAAGG")
    seq4 = LongSequence{DNAAlphabet{4}}("GAGGGGGATGCCCCTCTTTGAGCCCAAGG")
    
    #function test_iter_traits(it::CanonicalSkipmers{ST,U,S}) where {ST,U,S}
    function test_iter_traits(it::SkipmerFactory{S,U,M,N,K}) where {S,U,M,N,K}
        @test eltype(it) == SkipmerFactoryResult{U,DNAAlphabet{2},M,N,K}
        @test BioSequences.mertype(it) == Skipmer{U,DNAAlphabet{2},M,N,K}
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
    
    ST = Skipmer{UInt64, DNAAlphabet{2},2,3,14}
    BIGST = Skipmer{UInt128, DNAAlphabet{2},2,3,14}
    
    ans = collect(SkipmerFactoryResult(x...) for x in zip(eachindex(fwans), ST.(fwans), ST.(string_reverse_complement.(DNA, fwans))))
    bigans = collect(SkipmerFactoryResult(x...) for x in zip(eachindex(fwans), BIGST.(fwans), BIGST.(string_reverse_complement.(DNA, fwans))))
    
    #ans = ["CAAGGGGCTCCCTC", "AGGGATCCCTTTGA", "CTAAGAGGCACCCC",
    #       "GCCAAGGGGCTCCC", "GGATCCCTTTGACC", "GGCTAAGAGGCACC",
    #       "GAGCCCCTTGGCCA", "ATCCCTTTGACCAA", "CTGGCTAAGAGGCA",
    #       "CCTGGCCAAGGGGC"]

    #test_iter_traits(eachcanonical(ST, seq2))
    test_iter_traits(SkipmerFactory(ST, seq2))
    #test_iter_traits(eachcanonical(ST, seq4))
    test_iter_traits(SkipmerFactory(ST, seq4))
    #@test collect(eachcanonical(ST, seq2)) == ST.(ans)
    @test collect(SkipmerFactory(ST, seq2)) == ans
    #@test collect(eachcanonical(ST, seq4)) == ST.(ans)
    @test collect(SkipmerFactory(ST, seq4)) == ans
    
    # Test the same pattern, but with UInt128 for big capacity skipmers.
    #test_iter_traits(eachcanonical(BIGST, seq2))
    test_iter_traits(SkipmerFactory(BIGST, seq2))
    #test_iter_traits(eachcanonical(BIGST, seq4))
    test_iter_traits(SkipmerFactory(BIGST, seq4))
    #@test collect(eachcanonical(BIGST, seq2)) == BIGST.(ans)
    @test collect(SkipmerFactory(BIGST, seq2)) == bigans
    #@test collect(eachcanonical(BIGST, seq4)) == BIGST.(ans)
    @test collect(SkipmerFactory(BIGST, seq4)) == bigans
    
    # Test when Skipmers are based on UInt64, and sample nucleotides using
    # M = 1, N = 3, K = 7 pattern.
    ST2 = Skipmer{UInt64,DNAAlphabet{2},1,3,7}
    BIGST2 = Skipmer{UInt128,DNAAlphabet{2},1,3,7}
    
    fwans2 = ["GGGGCCT", "AGACCTG", "GGTCTTA", "GGGCCTG",
             "GACCTGC", "GTCTTAC", "GGCCTGC", "ACCTGCA",
             "TCTTACA", "GCCTGCG", "CCTGCAG"]
    
    ans2 = collect(zip(eachindex(fwans2), ST2.(fwans2), ST2.(string_reverse_complement.(DNA, fwans2))))
    bigans2 = collect(zip(eachindex(fwans2), BIGST2.(fwans2), BIGST2.(string_reverse_complement.(DNA, fwans2))))
    
    # Does iterator stop you if the span of your skipmer exceeds the len of the
    # sequence?
    @test_throws ArgumentError SkipmerFactory(Skipmer{UInt64, DNAAlphabet{2},1,3,14}, seq2)
    
    test_iter_traits(SkipmerFactory(ST2, seq2))
    test_iter_traits(SkipmerFactory(ST2, seq4))
    @test collect(SkipmerFactory(ST2, seq2)) == ans2
    @test collect(SkipmerFactory(ST2, seq4)) == ans2
    
    # Test the same pattern, but with UInt128 for big capacity skipmers.
    test_iter_traits(SkipmerFactory(BIGST2, seq2))
    test_iter_traits(SkipmerFactory(BIGST2, seq4))
    @test collect(SkipmerFactory(BIGST2, seq2)) == bigans2
    @test collect(SkipmerFactory(BIGST2, seq4)) == bigans2
    
    # Test when Skipmers are generated using M = 2, N = 4, K = 8 pattern.
    ST3 = Skipmer{UInt64, DNAAlphabet{2}, 2, 4, 8}
    BIGST3 = Skipmer{UInt128, DNAAlphabet{2}, 2, 4, 8}
    
    fwans3 = ["GAGGTGCC", "AGGGGCCT", "GGGACCTC", "GGATCCCT",
              "GGTGCCTT", "GGGCCTTT", "GACCTCTG", "ATCCCTGA",
              "TGCCTTAG", "GCCTTTGC", "CCTCTGCC","CCCTGACC",
              "CCTTAGCA", "CTTTGCAA", "TCTGCCAG", "CTGACCGG"]
    
    ans3 = collect(zip(eachindex(fwans3), ST3.(fwans3), ST3.(string_reverse_complement.(DNA, fwans3))))
    bigans3 = collect(zip(eachindex(fwans3), BIGST3.(fwans3), BIGST3.(string_reverse_complement.(DNA, fwans3))))
              
    test_iter_traits(SkipmerFactory(ST3, seq2))
    test_iter_traits(SkipmerFactory(ST3, seq4))
    @test collect(SkipmerFactory(ST3, seq2)) == ans3
    @test collect(SkipmerFactory(ST3, seq4)) == ans3
    
    # Test the same pattern, but with UInt128 for big capacity skipmers.
    test_iter_traits(SkipmerFactory(BIGST3, seq2))
    test_iter_traits(SkipmerFactory(BIGST3, seq4))
    @test collect(SkipmerFactory(BIGST3, seq2)) == bigans3
    @test collect(SkipmerFactory(BIGST3, seq4)) == bigans3
end

