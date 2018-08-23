@testset "Basic Operations" begin
    seq = ReferenceSequence(dna"")
    @test length(seq) === lastindex(seq) === 0
    @test isempty(seq)
    @test ReferenceSequence(dna"") == seq
    @test_throws BoundsError seq[0]
    @test_throws BoundsError seq[1]
    @test_throws BoundsError seq[2:3]
    @test collect(seq) == DNA[]

    seq = ReferenceSequence(dna"ACGTN")
    @test length(seq) === lastindex(seq) === 5
    @test !isempty(seq)
    @test seq[1] === DNA_A
    @test seq[2] === DNA_C
    @test seq[3] === DNA_G
    @test seq[4] === DNA_T
    @test seq[5] === DNA_N
    @test seq[2:3] == dna"CG"
    @test seq[3:5] == dna"GTN"
    @test seq[5:4] == dna""
    @test_throws BoundsError seq[0]
    @test_throws BoundsError seq[6]
    @test_throws BoundsError seq[0:2]
    @test_throws BoundsError seq[5:6]
    @test collect(seq) == [DNA_A, DNA_C, DNA_G, DNA_T, DNA_N]
end
