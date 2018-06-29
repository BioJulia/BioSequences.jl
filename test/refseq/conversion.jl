@testset "Construction" begin
    @test isa(ReferenceSequence(dna""), ReferenceSequence)
    @test isa(ReferenceSequence(dna"ACGTN"^30), ReferenceSequence)
    @test isa(ReferenceSequence("ACGTN"^30), ReferenceSequence)
    @test_throws Exception ReferenceSequence(dna"ACGRT")
    @test_throws Exception ReferenceSequence("ACGRT")
end

@testset "Conversion" begin
    seq = DNASequence(random_dna(100))
    refseq = ReferenceSequence(seq)
    @test refseq == seq
    @test DNASequence(refseq) == seq
end
