@testset "Print" begin
    buf = IOBuffer()
    print(buf, ReferenceSequence(dna""))
    @test String(take!(buf)) == ""

    print(buf, ReferenceSequence(dna"ACGTN"))
    @test String(take!(buf)) == "ACGTN"

    print(buf, ReferenceSequence(dna"A"^100))
    @test String(take!(buf)) == "A"^100

    print(buf, ReferenceSequence(dna"A"^100), width=80)
    @test String(take!(buf)) == "A"^80 * '\n' * "A"^20

    print(buf, ReferenceSequence(dna"A"^20), width=20)
    @test String(take!(buf)) == "A"^20

    # Check more than 4096
    print(buf, ReferenceSequence(dna"A"^4101), width=100)
    @test String(take!(buf)) == repeat("A"^100 * '\n', 41) * "A"
end

@testset "Show" begin
    buf = IOBuffer()
    show(buf, ReferenceSequence(dna""))
    @test String(take!(buf)) == "< EMPTY SEQUENCE >"

    show(buf, ReferenceSequence(dna"ACGTN"))
    @test String(take!(buf)) == "ACGTN"
end
