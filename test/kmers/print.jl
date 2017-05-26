@testset "Print" begin
    buf = IOBuffer()
    print(buf, DNAKmer(""))
    @test takebuf_string(buf) == ""

    buf = IOBuffer()
    print(buf, DNAKmer("ACGT"))
    @test takebuf_string(buf) == "ACGT"

    buf = IOBuffer()
    print(buf, RNAKmer("ACGU"))
    @test takebuf_string(buf) == "ACGU"
end
