@testset "Print" begin
    buf = IOBuffer()

    print(buf, DNAKmer("ACGT"))
    @test String(take!(buf)) == "ACGT"

    print(buf, RNAKmer("ACGU"))
    @test String(take!(buf)) == "ACGU"
end

@testset "Show" begin
    buf = IOBuffer()

    show(buf, DNAKmer("AGAGT"))
    @test String(take!(buf)) == "AGAGT"
end
