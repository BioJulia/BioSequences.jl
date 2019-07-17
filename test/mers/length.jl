@testset "Length" begin
    for len in [1, 16, 32]
        @test length(DNAMer(random_dna_kmer(len))) == len
        @test length(RNAMer(random_rna_kmer(len))) == len
    end
end
