@testset "Shuffle" begin
    @testset for i in 1:10
        @test shuffle(dna"") == dna""
        @test shuffle(dna"A") == dna"A"
        @test shuffle(dna"C") == dna"C"
    end

    seq = dna"ACGTN"^10
    @test shuffle(seq) != dna"ACGTN"^10
    @test seq == dna"ACGTN"^10
    @test shuffle!(seq) === seq
    @test seq != dna"ACGTN"^10
    @test count(x -> x == DNA_A, seq) == 10
    @test count(x -> x == DNA_C, seq) == 10
    @test count(x -> x == DNA_G, seq) == 10
    @test count(x -> x == DNA_T, seq) == 10
    @test count(x -> x == DNA_N, seq) == 10
end
