@testset "MatchQuery" begin
    seq = dna"ACGTACG"
    
    @testset "forward" begin
        @test findfirst(MatchQuery(dna""), seq) === 1:0
        @test findfirst(MatchQuery(dna"AC"), seq) === 1:2
        @test findnext(MatchQuery(dna"AC"), seq, 2) === 5:6
        #@test findfirst(MatchQuery(dna"AC"), seq, 2, 5) === nothing
        @test findfirst(MatchQuery(dna"TG"), seq) === nothing
        @test findfirst(MatchQuery(dna"TN"), seq) === nothing
        @test findfirst(MatchQuery(dna"ACG"), seq)  === 1:3
        @test findnext(MatchQuery(dna"ACG"), seq, 2) === 5:7
        @test findfirst(MatchQuery(seq), seq) === 1:lastindex(seq)

        @test findfirst(MatchQuery(dna""), dna"") === 1:0
        @test findnext(MatchQuery(dna""), dna"", -1) === 1:0
        @test findnext(MatchQuery(dna""), dna"", 2) === nothing

        @test first(findfirst(MatchQuery(dna""), seq)) === 1
        @test first(findfirst(MatchQuery(dna"AC"), seq)) === 1
        @test first(findnext(MatchQuery(dna"AC"), seq, 2)) === 5
        #@test findnext(MatchQuery(dna"AC"), seq, 2) === nothing

        query = MatchQuery(dna"ACG")
        @test findfirst(query, seq) === 1:3
        @test findnext(query, seq, 2) === 5:7
        #@test findfirst(query, seq, 2, 6) === nothing
        @test first(findfirst(query, seq)) === 1
        @test first(findnext(query, seq, 2)) === 5
        #@test findfirst(query, seq, 2, 6) === nothing
    end

    @testset "backward" begin
        @test findlast(MatchQuery(dna""), seq) === 8:7
        @test findlast(MatchQuery(dna"AC"), seq) === 5:6
        @test findprev(MatchQuery(dna"AC"), seq, 5) === 1:2
        #@test findlast(MatchQuery(dna"AC"), seq, 5, 2) === nothing
        @test findlast(MatchQuery(dna"TG"), seq) === nothing
        @test findlast(MatchQuery(dna"TN"), seq) === nothing
        @test findlast(MatchQuery(dna"ACG"), seq) === 5:7
        @test findprev(MatchQuery(dna"ACG"), seq, 6) === 1:3
        @test findlast(MatchQuery(seq), seq) === 1:lastindex(seq)

        @test findlast(MatchQuery(dna""), dna"") === 1:0
        @test findprev(MatchQuery(dna""), dna"", 2) === 1:0
        @test findprev(MatchQuery(dna""), dna"", -1) === nothing

        @test first(findlast(MatchQuery(dna""), seq)) === 8
        @test first(findlast(MatchQuery(dna"AC"), seq)) === 5
        @test first(findprev(MatchQuery(dna"AC"), seq, 5)) === 1
        #@test findlast(MatchQuery(dna"AC"), seq, 5, 2) === nothing

        query = MatchQuery(dna"ACG")
        @test findlast(query, seq) === 5:7
        @test findprev(query, seq, 6) === 1:3
        #@test findlast(query, seq, 2, 6) === nothing
        @test first(findlast(query, seq)) === 5
        @test first(findprev(query, seq, 6)) === 1
        #@test findlast(query, seq, 6, 2) === nothing
    end

    @testset "occursin" begin
        @test occursin(MatchQuery(dna"ACG"), dna"GGGTACACGTTT") == true
        @test occursin(MatchQuery(dna"TGT"), dna"GGGTACACGTGT") == true
        @test occursin(MatchQuery(dna"GGG"), dna"GGGTACACGTGT") == true
        @test occursin(MatchQuery(dna"AAA"), dna"GGGTACACGTGT") == false
    end

end
