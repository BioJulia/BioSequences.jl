@testset "Each k-mer" begin
    function string_eachkmer(seq::AbstractString, k, step)
        kmers = String[]
        i = 1
        for i in 1:step:length(seq) - k + 1
            subseq = seq[i:i + k - 1]
            if !in('N', subseq)
                push!(kmers, subseq)
            end
        end
        return kmers
    end

    function test_eachkmer(S, seq::AbstractString, k, step)
        xs = [convert(String, x)
              for (i, x) in collect(each(Kmer{eltype(S),k}, S(seq), step))]
        ys = [convert(String, x)
              for (i, x) in collect(eachkmer(S(seq), k, step))]
        zs = string_eachkmer(seq, k, step)
        @test xs == ys == zs
    end

    for k in [0, 1, 3, 16, 32], step in 1:3, len in [1, 2, 3, 5, 10, 100, 1000]
        test_eachkmer(BioSequence{DNAAlphabet{4}}, random_dna(len), k, step)
        test_eachkmer(BioSequence{RNAAlphabet{4}}, random_rna(len), k, step)
        test_eachkmer(ReferenceSequence, random_dna(len), k, step)

        probs = [0.25, 0.25, 0.25, 0.25, 0.00]
        test_eachkmer(BioSequence{DNAAlphabet{2}}, random_dna(len, probs), k, step)
        test_eachkmer(BioSequence{RNAAlphabet{2}}, random_rna(len, probs), k, step)
    end

    @test isempty(collect(each(DNAKmer{1}, dna"")))
    @test isempty(collect(each(DNAKmer{1}, dna"NNNNNNNNNN")))
    @test_throws Exception each(DNAKmer{-1}, dna"ACGT")
    @test_throws Exception each(DNAKmer{33}, dna"ACGT")
end
