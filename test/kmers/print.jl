@testset "Print" begin
    buf = IOBuffer()
    print(buf, DNAKmer(""))
    @test String(take!(buf)) == ""

    buf = IOBuffer()
    print(buf, DNAKmer("ACGT"))
    @test String(take!(buf)) == "ACGT"

    buf = IOBuffer()
    print(buf, RNAKmer("ACGU"))
    @test String(take!(buf)) == "ACGU"
end
