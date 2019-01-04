@testset "Transformations" begin
    function test_reverse(T, seq)
        revseq = reverse(Kmer{T,length(seq)}(seq))
        @test String(revseq) == reverse(seq)
    end

    function test_dna_complement(seq)
        comp = complement(DNAKmer{length(seq)}(seq))
        @test String(comp) == dna_complement(seq)
    end

    function test_rna_complement(seq)
        comp = complement(RNAKmer{length(seq)}(seq))
        @test String(comp) == rna_complement(seq)
    end

    function test_dna_revcomp(seq)
        revcomp = reverse_complement(DNAKmer{length(seq)}(seq))
        @test String(revcomp) == reverse(dna_complement(seq))
    end

    function test_rna_revcomp(seq)
        revcomp = reverse_complement(RNAKmer{length(seq)}(seq))
        @test String(revcomp) == reverse(rna_complement(seq))
    end

    @testset "Reverse" begin
        for len in 1:32, _ in 1:10
            test_reverse(DNAAlphabet{2}, random_dna_kmer(len))
            test_reverse(RNAAlphabet{2}, random_rna_kmer(len))
        end

        seq = dna"AAAAAAAAAAAAAAAAAAAAAAAAAAAAGATAC"
        @test reverse(seq[(length(seq)-9):length(seq)]) == dna"CATAGAAAAA"
    end

    @testset "Complement" begin
        for len in 1:32, _ in 1:10
            test_dna_complement(random_dna_kmer(len))
            test_rna_complement(random_rna_kmer(len))
        end
    end

    @testset "Reverse Complement" begin
        for len in 1:32, _ in 1:10
            test_dna_revcomp(random_dna_kmer(len))
            test_rna_revcomp(random_rna_kmer(len))
        end
    end
end
