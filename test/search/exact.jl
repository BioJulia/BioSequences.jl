@testset "SearchQuery" begin
    seq = dna"ACGTACG"
    
    @testset "forward" begin
        @test findfirst(SearchQuery(dna""), seq) === 1:0
        @test findfirst(SearchQuery(dna"AC"), seq) === 1:2
        @test findnext(SearchQuery(dna"AC"), seq, 2) === 5:6
        #@test findfirst(SearchQuery(dna"AC"), seq, 2, 5) === nothing
        @test findfirst(SearchQuery(dna"TG"), seq) === nothing
        @test findfirst(SearchQuery(dna"TN"), seq) === 4:5
        @test findfirst(SearchQuery(dna"ACG"), seq)  === 1:3
        @test findnext(SearchQuery(dna"ACG"), seq, 2) === 5:7
        @test findfirst(SearchQuery(seq), seq) === 1:lastindex(seq)

        @test findfirst(DNA_A, dna"AAACT") == 1
        @test findfirst(DNA_C, dna"") === nothing
        @test findfirst(DNA_M, dna"TGC") == 3
        @test findfirst(DNA_N, dna"T") == 1
        @test findfirst(DNA_V, dna"TTTT") === nothing
        @test findfirst(DNA_C, seq) == 2

        @test findfirst(dna"", dna"") === 1:0
        @test findfirst(dna"", dna"", -1) === 1:0
        @test findfirst(dna"", dna"", 2) === nothing

        @test first(findfirst(SearchQuery(dna""), seq)) === 1
        @test first(findfirst(SearchQuery(dna"AC"), seq)) === 1
        @test first(findnext(SearchQuery(dna"AC"), seq, 2)) === 5
        #@test findnext(SearchQuery(dna"AC"), seq, 2) === nothing

        query = SearchQuery(dna"ACG")
        @test findfirst(query, seq) === 1:3
        @test findnext(query, seq, 2) === 5:7
        #@test findfirst(query, seq, 2, 6) === nothing
        @test first(findfirst(query, seq)) === 1
        @test first(findnext(query, seq, 2)) === 5
        #@test findfirst(query, seq, 2, 6) === nothing
    end

    @testset "backward" begin
        @test findlast(SearchQuery(dna""), seq) === 8:7
        @test findlast(SearchQuery(dna"AC"), seq) === 5:6
        @test findprev(SearchQuery(dna"AC"), seq, 5) === 1:2
        #@test findlast(SearchQuery(dna"AC"), seq, 5, 2) === nothing
        @test findlast(SearchQuery(dna"TG"), seq) === nothing
        @test findlast(SearchQuery(dna"TN"), seq) === 4:5
        @test findlast(SearchQuery(dna"ACG"), seq) === 5:7
        @test findprev(SearchQuery(dna"ACG"), seq, 6) === 1:3
        @test findlast(SearchQuery(seq), seq) === 1:lastindex(seq)

        @test findlast(DNA_A, seq) == 5
        @test findlast(DNA_T, seq) == 4
        @test findlast(DNA_N, dna"") === nothing
        @test findlast(DNA_R, dna"AGTCAGTTCT") == 6
        @test findlast(DNA_R, dna"TCTCTT") === nothing

        @test findlast(dna"", dna"") === 1:0
        @test findlast(dna"", dna"", 2) === 1:0
        @test findlast(dna"", dna"", -1) === nothing

        @test first(findlast(SearchQuery(dna""), seq)) === 8
        @test first(findlast(SearchQuery(dna"AC"), seq)) === 5
        @test first(findprev(SearchQuery(dna"AC"), seq, 5)) === 1
        #@test findlast(SearchQuery(dna"AC"), seq, 5, 2) === nothing

        query = SearchQuery(dna"ACG")
        @test findlast(query, seq) === 5:7
        @test findprev(query, seq, 6) === 1:3
        #@test findlast(query, seq, 2, 6) === nothing
        @test first(findlast(query, seq)) === 5
        @test first(findprev(query, seq, 6)) === 1
        #@test findlast(query, seq, 6, 2) === nothing
    end

    @testset "occursin" begin
        @test occursin(SearchQuery(dna"ACG"), dna"GGGTACACGTTT") == true
        @test occursin(SearchQuery(dna"TGT"), dna"GGGTACACGTGT") == true
        @test occursin(SearchQuery(dna"GGG"), dna"GGGTACACGTGT") == true
        @test occursin(SearchQuery(dna"AAA"), dna"GGGTACACGTGT") == false
    end

end
