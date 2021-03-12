@testset "Transformations" begin
    function test_reverse(T, seq)
        revseq = reverse(T(seq))
        @test String(revseq) == reverse(seq)
    end

    function test_dna_complement(T, seq)
        comp = complement(T(seq))
        @test String(comp) == dna_complement(seq)
    end

    function test_rna_complement(T, seq)
        comp = complement(T(seq))
        @test String(comp) == rna_complement(seq)
    end

    function test_dna_revcomp(T, seq)
        revcomp = reverse_complement(T(seq))
        @test String(revcomp) == reverse(dna_complement(seq))
    end

    function test_rna_revcomp(T, seq)
        revcomp = reverse_complement(T(seq))
        @test String(revcomp) == reverse(rna_complement(seq))
    end

    @testset "Reverse" begin
        for len in [1, 3, 5, 8, 15, 16, 30, 32, 55, 64], _ in 1:5
            if len <= 32
                test_reverse(DNAMer{len}, random_dna_kmer(len))
                test_reverse(RNAMer{len}, random_rna_kmer(len))
            end
            test_reverse(BigDNAMer{len}, random_dna_kmer(len))
            test_reverse(BigRNAMer{len}, random_rna_kmer(len))
        end

        seq = dna"AAAAAAAAAAAAAAAAAAAAAAAAAAAAGATAC"
        @test reverse(seq[(length(seq)-9):length(seq)]) == dna"CATAGAAAAA"
    end

    @testset "Complement" begin
        for len in [1, 3, 5, 8, 15, 16, 30, 32, 55, 64], _ in 1:5
            if len <= 32
                test_dna_complement(DNAMer{len}, random_dna_kmer(len))
                test_rna_complement(RNAMer{len}, random_rna_kmer(len))
            end
            test_dna_complement(BigDNAMer{len}, random_dna_kmer(len))
            test_rna_complement(BigRNAMer{len}, random_rna_kmer(len))
        end
    end

    @testset "Reverse Complement" begin
        for len in [1, 3, 5, 8, 15, 16, 30, 32, 55, 64], _ in 1:5
            if len <= 32
                test_dna_revcomp(DNAMer{len}, random_dna_kmer(len))
                test_rna_revcomp(RNAMer{len}, random_rna_kmer(len))
            end
            test_dna_revcomp(BigDNAMer{len}, random_dna_kmer(len))
            test_rna_revcomp(BigRNAMer{len}, random_rna_kmer(len))
        end
    end
end
