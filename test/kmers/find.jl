@testset "Find" begin
    kmer = DNAKmer("ACGAG")

    @test findnext(kmer, DNA_A, 1) == 1
    @test findnext(kmer, DNA_C, 1) == 2
    @test findnext(kmer, DNA_G, 1) == 3
    @test findnext(kmer, DNA_T, 1) == 0
    @test findnext(kmer, DNA_A, 2) == 4

    @test_throws BoundsError findnext(kmer, DNA_A, 0)
    @test_throws BoundsError findnext(kmer, DNA_A, 6)

    @test findprev(kmer, DNA_A, 5) == 4
    @test findprev(kmer, DNA_C, 5) == 2
    @test findprev(kmer, DNA_G, 5) == 5
    @test findprev(kmer, DNA_T, 5) == 0
    @test findprev(kmer, DNA_G, 4) == 3

    @test_throws BoundsError findprev(kmer, DNA_A, 0)
    @test_throws BoundsError findprev(kmer, DNA_A, 6)

    @test findfirst(kmer, DNA_A) == 1
    @test findfirst(kmer, DNA_G) == 3
    @test findlast(kmer, DNA_A) == 4
    @test findlast(kmer, DNA_G) == 5
end
