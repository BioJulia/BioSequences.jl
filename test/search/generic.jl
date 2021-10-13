@testset "findall" begin

    seq = dna"NNNN" * dna"GATC" * dna"GATC"

    r = findall(dna"GATC", seq)
    @test r == [1:4, 5:8, 9:12]

    r = findall(isequal(dna"GATC"), seq)
    @test r == [5:8, 9:12]

    # Check start advancement.
    @test [findfirst(dna"A", dna"GATC")] == findall(dna"A", dna"GATC")

    # Check start/stop delimiter.
    @test findall(dna"A", seq * dna"AA", r) == findall(dna"A", seq)

    # Check vcat of results.
    @test findall(dna"GATC", seq, r[1]) == findall(dna"GATC", seq, [r[1]])

end # testset findall
