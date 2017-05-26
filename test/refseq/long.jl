@testset "Long Sequence" begin
    str = random_dna(50000)
    @test ReferenceSequence(str) == DNASequence(str)
    str = ("N"^300 * random_dna(1000))^5 * "N"^300
    @test ReferenceSequence(str) == DNASequence(str)
end
