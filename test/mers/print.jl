@testset "Print" begin
    buf = IOBuffer()

    print(buf, DNAMer("ACGT"))
    @test String(take!(buf)) == "ACGT"

    print(buf, RNAMer("ACGU"))
    @test String(take!(buf)) == "ACGU"
end

@testset "Show" begin
    buf = IOBuffer()

    show(buf, DNAMer("AGAGT"))
    @test String(take!(buf)) == "AGAGT"
end
