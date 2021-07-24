@testset "Exact" begin
    seq = dna"ACGTACG"

    @testset "forward" begin
        @test findfirst(dna"", seq) === 1:0
        @test findfirst(dna"AC", seq) === 1:2
        @test findfirst(dna"AC", seq, 2) === 5:6
        @test findfirst(dna"AC", seq, 2, 5) === nothing
        @test findfirst(dna"TG", seq) === nothing
        @test findfirst(dna"TN", seq) === 4:5
        @test findfirst(dna"ACG", seq)  === 1:3
        @test findfirst(dna"ACG", seq, 2) === 5:7
        @test findfirst(seq, seq) === 1:lastindex(seq)

        @test findfirst(DNA_A, dna"AAACT") == 1
        @test findfirst(DNA_C, dna"") === nothing
        @test findfirst(DNA_M, dna"TGC") == 3
        @test findfirst(DNA_N, dna"T") == 1
        @test findfirst(DNA_V, dna"TTTT") === nothing
        @test findfirst(DNA_C, seq) == 2

        @test findfirst(dna"", dna"") === 1:0
        @test findfirst(dna"", dna"", -1) === 1:0
        @test findfirst(dna"", dna"", 2) === nothing

        @test first(findfirst(dna"", seq)) === 1
        @test first(findfirst(dna"AC", seq)) === 1
        @test first(findfirst(dna"AC", seq, 2)) === 5
        @test findfirst(dna"AC", seq, 2, 5) === nothing

        query = ExactSearchQuery(dna"ACG")
        @test findfirst(query, seq) === 1:3
        @test findfirst(query, seq, 2) === 5:7
        @test findfirst(query, seq, 2, 6) === nothing
        @test first(findfirst(query, seq)) === 1
        @test first(findfirst(query, seq, 2)) === 5
        @test findfirst(query, seq, 2, 6) === nothing
    end

    @testset "backward" begin
        @test findlast(dna"", seq) === 8:7
        @test findlast(dna"AC", seq) === 5:6
        @test findlast(dna"AC", seq, 5) === 1:2
        @test findlast(dna"AC", seq, 5, 2) === nothing
        @test findlast(dna"TG", seq) === nothing
        @test findlast(dna"TN", seq) === 4:5
        @test findlast(dna"ACG", seq) === 5:7
        @test findlast(dna"ACG", seq, 6) === 1:3
        @test findlast(seq, seq) === 1:lastindex(seq)

        @test findlast(DNA_A, seq) == 5
        @test findlast(DNA_T, seq) == 4
        @test findlast(DNA_N, dna"") === nothing
        @test findlast(DNA_R, dna"AGTCAGTTCT") == 6
        @test findlast(DNA_R, dna"TCTCTT") === nothing

        @test findlast(dna"", dna"") === 1:0
        @test findlast(dna"", dna"", 2) === 1:0
        @test findlast(dna"", dna"", -1) === nothing

        @test first(findlast(dna"", seq)) === 8
        @test first(findlast(dna"AC", seq)) === 5
        @test first(findlast(dna"AC", seq, 5)) === 1
        @test findlast(dna"AC", seq, 5, 2) === nothing

        query = ExactSearchQuery(dna"ACG")
        @test findlast(query, seq) === 5:7
        @test findlast(query, seq, 6) === 1:3
        @test findlast(query, seq, 2, 6) === nothing
        @test first(findlast(query, seq)) === 5
        @test first(findlast(query, seq, 6)) === 1
        @test findlast(query, seq, 6, 2) === nothing
    end

    @testset "occursin" begin
        @test occursin(dna"ACG", dna"GGGTACACGTTT") == true
        @test occursin(dna"TGT", dna"GGGTACACGTGT") == true
        @test occursin(dna"GGG", dna"GGGTACACGTGT") == true
        @test occursin(dna"AAA", dna"GGGTACACGTGT") == false
    end

end
