@testset "Length" begin
    for len in [1, 16, 32]
        @test length(DNAKmer(random_dna_kmer(len))) == len
        @test length(RNAKmer(random_rna_kmer(len))) == len
    end
end
