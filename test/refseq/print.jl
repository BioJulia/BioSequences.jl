@testset "Print" begin
    buf = IOBuffer()
    print(buf, ReferenceSequence(dna""))
    @test String(take!(buf)) == ""

    print(buf, ReferenceSequence(dna"ACGTN"))
    @test String(take!(buf)) == "ACGTN"

    print(buf, ReferenceSequence(dna"A"^100))
    @test String(take!(buf)) == "A"^100
end

@testset "Show" begin
    buf = IOBuffer()
    show(buf, ReferenceSequence(dna""))
    @test String(take!(buf)) == "< EMPTY SEQUENCE >"

    show(buf, ReferenceSequence(dna"ACGTN"))
    @test String(take!(buf)) == "ACGTN"
end
