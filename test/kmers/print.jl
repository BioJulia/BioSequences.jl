@testset "Print" begin
    buf = IOBuffer()

    print(buf, DNAKmer(""))
    @test String(take!(buf)) == ""

    print(buf, DNAKmer("ACGT"))
    @test String(take!(buf)) == "ACGT"

    print(buf, RNAKmer("ACGU"))
    @test String(take!(buf)) == "ACGU"
end

@testset "Show" begin
    buf = IOBuffer()

    show(buf, DNAKmer(""))
    @test String(take!(buf)) == "DNA 0-mer:\n< EMPTY SEQUENCE >"

    show(buf, DNAKmer("AGAGT"))
    @test String(take!(buf)) == "DNA 5-mer:\nAGAGT"
end
