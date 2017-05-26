@testset "Find" begin
    seq = dna"ACGNA"
    @test findnext(seq, DNA_A, 1) == 1
    @test findnext(seq, DNA_C, 1) == 2
    @test findnext(seq, DNA_G, 1) == 3
    @test findnext(seq, DNA_N, 1) == 4
    @test findnext(seq, DNA_T, 1) == 0
    @test findnext(seq, DNA_A, 2) == 5

    @test_throws BoundsError findnext(seq, DNA_A, 0)
    @test_throws BoundsError findnext(seq, DNA_A, 6)

    @test findprev(seq, DNA_A, 4) == 1
    @test findprev(seq, DNA_C, 4) == 2
    @test findprev(seq, DNA_G, 4) == 3
    @test findprev(seq, DNA_N, 4) == 4
    @test findprev(seq, DNA_T, 4) == 0
    @test findprev(seq, DNA_G, 2) == 0

    @test_throws BoundsError findprev(seq, DNA_A, 0)
    @test_throws BoundsError findprev(seq, DNA_A, 6)

    seq = dna"ACGNAN"
    @test findfirst(seq, DNA_A) == 1
    @test findfirst(seq, DNA_N) == 4
    @test findfirst(seq, DNA_T) == 0

    @test findlast(seq, DNA_A) == 5
    @test findlast(seq, DNA_N) == 6
    @test findlast(seq, DNA_T) == 0
end
