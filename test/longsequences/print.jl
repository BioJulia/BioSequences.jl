@testset "Print" begin
    buf = IOBuffer()

    print(buf, dna"")
    @test String(take!(buf)) == ""

    print(buf, dna"ACGTN")
    @test String(take!(buf)) == "ACGTN"

    print(buf, rna"ACGUN")
    @test String(take!(buf)) == "ACGUN"

    print(buf, dna"A"^100)
    @test String(take!(buf)) == "A"^100

    print(buf, dna"G"^60, width=0)
    @test String(take!(buf)) == "G"^60

    print(buf, dna"A"^60, width=-5)
    @test String(take!(buf)) == "A"^60

    print(buf, dna"A"^100, width=70)
    @test String(take!(buf)) == string("A"^70, '\n', "A"^30)

    print(buf, dna"A"^100, width=50)
    @test String(take!(buf)) == string("A"^50, '\n', "A"^50)

    print(buf, dna"A"^4100, width=100)
    @test String(take!(buf)) == repeat("A"^100 * '\n', 40) * "A"^100
end

@testset "Show" begin
    buf = IOBuffer()

    show(buf, dna"")
    @test String(take!(buf)) == "< EMPTY SEQUENCE >"

    show(buf, dna"ATCG")
    @test String(take!(buf)) == "ATCG"
end
