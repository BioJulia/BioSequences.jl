@testset "Find" begin
    seq = dna"ACGNA"
    @test findnext(isequal(DNA_A), seq, 1) == 1
    @test findnext(isequal(DNA_C), seq, 1) == 2
    @test findnext(isequal(DNA_G), seq, 1) == 3
    @test findnext(isequal(DNA_N), seq, 1) == 4
    @test findnext(isequal(DNA_T), seq, 1) === nothing
    @test findnext(isequal(DNA_A), seq, 2) == 5

    @test_throws BoundsError findnext(isequal(DNA_A), seq, 0)
    @test findnext(isequal(DNA_A), seq, 6) === nothing

    @test findprev(isequal(DNA_A), seq, 4) == 1
    @test findprev(isequal(DNA_C), seq, 4) == 2
    @test findprev(isequal(DNA_G), seq, 4) == 3
    @test findprev(isequal(DNA_N), seq, 4) == 4
    @test findprev(isequal(DNA_T), seq, 4) === nothing
    @test findprev(isequal(DNA_G), seq, 2) === nothing

    @test findprev(isequal(DNA_A), seq, 0) === nothing
    @test_throws BoundsError findprev(isequal(DNA_A), seq, 6)

    seq = dna"ACGNAN"
    @test findfirst(isequal(DNA_A), seq) == 1
    @test findfirst(isequal(DNA_N), seq) == 4
    @test findfirst(isequal(DNA_T), seq) === nothing

    @test findlast(isequal(DNA_A), seq) == 5
    @test findlast(isequal(DNA_N), seq) == 6
    @test findlast(isequal(DNA_T), seq) === nothing


    #         0000000001111
    #         1234567890123
    seq = dna"NNNNNGATCGATC"

    # Check default search.
    @test findall(DNA_A, seq) == [7, 11]
    @test findall(ExactSearchQuery(dna"A"), seq) == [7:7, 11:11]

    # Check overlap key argument.
    @test findall(ExactSearchQuery(dna"GATC", iscompatible), seq; overlap = false) == [1:4, 6:9, 10:13]
    @test findall(ExactSearchQuery(dna"GATC", iscompatible), seq; overlap = true) == [1:4, 2:5, 6:9, 10:13]

    # Check mapping of indices.
    @test findall(DNA_A, seq, 7:11) == [7, 11]
    @test findall(ExactSearchQuery(dna"A"), seq, 7:11) == [7:7, 11:11]

    # Check empty return type.
    @test findall(DNA_A, dna"GGGG") |> typeof == Vector{Int}
    @test findall(ExactSearchQuery(dna"A"), dna"GGGG") |> typeof == Vector{UnitRange{Int}}

end
