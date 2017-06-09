@testset "Print" begin
    buf = IOBuffer()
    print(buf, dna"")
    @test String(take!(buf)) == ""

    buf = IOBuffer()
    print(buf, dna"ACGTN")
    @test String(take!(buf)) == "ACGTN"

    buf = IOBuffer()
    print(buf, rna"ACGUN")
    @test String(take!(buf)) == "ACGUN"

    buf = IOBuffer()
    print(buf, dna"A"^100)
    @test String(take!(buf)) == "A"^100

    buf = IOBuffer()
    print(buf, dna"A"^100, width=70)
    @test String(take!(buf)) == string("A"^70, '\n', "A"^30)
end
