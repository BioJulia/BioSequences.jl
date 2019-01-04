@testset "Find" begin
    seq = dna"ACGNA"
    @test findnext(DNA_A, seq, 1) == 1
    @test findnext(DNA_C, seq, 1) == 2
    @test findnext(DNA_G, seq, 1) == 3
    @test findnext(DNA_N, seq, 1) == 4
    @test findnext(DNA_T, seq, 1) == nothing
    @test findnext(DNA_A, seq, 2) == 5

    @test_throws BoundsError findnext(DNA_A, seq, 0)
    @test_throws BoundsError findnext(DNA_A, seq, 6)

    @test findprev(DNA_A, seq, 4) == 1
    @test findprev(DNA_C, seq, 4) == 2
    @test findprev(DNA_G, seq, 4) == 3
    @test findprev(DNA_N, seq, 4) == 4
    @test findprev(DNA_T, seq, 4) == nothing
    @test findprev(DNA_G, seq, 2) == nothing

    @test_throws BoundsError findprev(DNA_A, seq, 0)
    @test_throws BoundsError findprev(DNA_A, seq, 6)

    seq = dna"ACGNAN"
    @test findfirst(DNA_A, seq) == 1
    @test findfirst(DNA_N, seq) == 4
    @test findfirst(DNA_T, seq) == nothing

    @test findlast(DNA_A, seq) == 5
    @test findlast(DNA_N, seq) == 6
    @test findlast(DNA_T, seq) == nothing
end
