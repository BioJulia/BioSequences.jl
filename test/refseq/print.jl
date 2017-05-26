@testset "Print" begin
    buf = IOBuffer()
    print(buf, ReferenceSequence(dna""))
    @test takebuf_string(buf) == ""

    buf = IOBuffer()
    print(buf, ReferenceSequence(dna"ACGTN"))
    @test takebuf_string(buf) == "ACGTN"

    buf = IOBuffer()
    print(buf, ReferenceSequence(dna"A"^100))
    @test takebuf_string(buf) == "A"^100
end
