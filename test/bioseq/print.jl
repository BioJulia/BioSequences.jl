@testset "Print" begin
    buf = IOBuffer()
    print(buf, dna"")
    @test takebuf_string(buf) == ""

    buf = IOBuffer()
    print(buf, dna"ACGTN")
    @test takebuf_string(buf) == "ACGTN"

    buf = IOBuffer()
    print(buf, rna"ACGUN")
    @test takebuf_string(buf) == "ACGUN"

    buf = IOBuffer()
    print(buf, dna"A"^100)
    @test takebuf_string(buf) == "A"^100

    buf = IOBuffer()
    print(buf, dna"A"^100, width=70)
    @test takebuf_string(buf) == string("A"^70, '\n', "A"^30)
end
