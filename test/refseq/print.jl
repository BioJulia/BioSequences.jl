@testset "Print" begin
    buf = IOBuffer()
    print(buf, ReferenceSequence(dna""))
    @test String(take!(buf)) == ""

    buf = IOBuffer()
    print(buf, ReferenceSequence(dna"ACGTN"))
    @test String(take!(buf)) == "ACGTN"

    buf = IOBuffer()
    print(buf, ReferenceSequence(dna"A"^100))
    @test String(take!(buf)) == "A"^100
end
