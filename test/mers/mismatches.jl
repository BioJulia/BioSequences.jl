@testset "Mismatches" begin
    function test_mismatches(a, b)
        count = 0
        for (x, y) in zip(a, b)
            count += x != y
        end
        @test mismatches(a, b) === mismatches(b, a) === count
    end

    for len in 1:64, _ in 1:10
        a = random_dna_kmer(len)
        b = random_dna_kmer(len)
        if len <= 32
            test_mismatches(DNAMer(a), DNAMer(b))
        end
        test_mismatches(BigDNAMer(a), BigDNAMer(b))
        
        a = random_rna_kmer(len)
        b = random_rna_kmer(len)
        if len <= 32
            test_mismatches(RNAMer(a), RNAMer(b))
        end
        test_mismatches(BigRNAMer(a), BigRNAMer(b))
    end
end
