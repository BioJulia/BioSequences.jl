@testset "MatchSearchQuery" begin
    seq = dna"ACGTACG"
    
    @testset "forward" begin
        @test findfirst(MatchSearchQuery(dna""), seq) === 1:0
        @test findfirst(MatchSearchQuery(dna"AC"), seq) === 1:2
        @test findnext(MatchSearchQuery(dna"AC"), seq, 2) === 5:6
        #@test findfirst(MatchSearchQuery(dna"AC"), seq, 2, 5) === nothing
        @test findfirst(MatchSearchQuery(dna"TG"), seq) === nothing
        @test findfirst(MatchSearchQuery(dna"TN"), seq) === nothing
        @test findfirst(MatchSearchQuery(dna"ACG"), seq)  === 1:3
        @test findnext(MatchSearchQuery(dna"ACG"), seq, 2) === 5:7
        @test findfirst(MatchSearchQuery(seq), seq) === 1:lastindex(seq)

        @test findfirst(MatchSearchQuery(dna""), dna"") === 1:0
        @test findnext(MatchSearchQuery(dna""), dna"", -1) === 1:0
        @test findnext(MatchSearchQuery(dna""), dna"", 2) === nothing

        @test first(findfirst(MatchSearchQuery(dna""), seq)) === 1
        @test first(findfirst(MatchSearchQuery(dna"AC"), seq)) === 1
        @test first(findnext(MatchSearchQuery(dna"AC"), seq, 2)) === 5
        #@test findnext(MatchSearchQuery(dna"AC"), seq, 2) === nothing

        query = MatchSearchQuery(dna"ACG")
        @test findfirst(query, seq) === 1:3
        @test findnext(query, seq, 2) === 5:7
        #@test findfirst(query, seq, 2, 6) === nothing
        @test first(findfirst(query, seq)) === 1
        @test first(findnext(query, seq, 2)) === 5
        #@test findfirst(query, seq, 2, 6) === nothing
    end

    @testset "backward" begin
        @test findlast(MatchSearchQuery(dna""), seq) === 8:7
        @test findlast(MatchSearchQuery(dna"AC"), seq) === 5:6
        @test findprev(MatchSearchQuery(dna"AC"), seq, 5) === 1:2
        #@test findlast(MatchSearchQuery(dna"AC"), seq, 5, 2) === nothing
        @test findlast(MatchSearchQuery(dna"TG"), seq) === nothing
        @test findlast(MatchSearchQuery(dna"TN"), seq) === nothing
        @test findlast(MatchSearchQuery(dna"ACG"), seq) === 5:7
        @test findprev(MatchSearchQuery(dna"ACG"), seq, 6) === 1:3
        @test findlast(MatchSearchQuery(seq), seq) === 1:lastindex(seq)

        @test findlast(MatchSearchQuery(dna""), dna"") === 1:0
        @test findprev(MatchSearchQuery(dna""), dna"", 2) === 1:0
        @test findprev(MatchSearchQuery(dna""), dna"", -1) === nothing

        @test first(findlast(MatchSearchQuery(dna""), seq)) === 8
        @test first(findlast(MatchSearchQuery(dna"AC"), seq)) === 5
        @test first(findprev(MatchSearchQuery(dna"AC"), seq, 5)) === 1
        #@test findlast(MatchSearchQuery(dna"AC"), seq, 5, 2) === nothing

        query = MatchSearchQuery(dna"ACG")
        @test findlast(query, seq) === 5:7
        @test findprev(query, seq, 6) === 1:3
        #@test findlast(query, seq, 2, 6) === nothing
        @test first(findlast(query, seq)) === 5
        @test first(findprev(query, seq, 6)) === 1
        #@test findlast(query, seq, 6, 2) === nothing
    end

    @testset "occursin" begin
        @test occursin(MatchSearchQuery(dna"ACG"), dna"GGGTACACGTTT") == true
        @test occursin(MatchSearchQuery(dna"TGT"), dna"GGGTACACGTGT") == true
        @test occursin(MatchSearchQuery(dna"GGG"), dna"GGGTACACGTGT") == true
        @test occursin(MatchSearchQuery(dna"AAA"), dna"GGGTACACGTGT") == false
    end

end
