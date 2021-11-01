@testset "CompatibleSearchQuery" begin
    seq = dna"ACGTACG"
    
    @testset "forward" begin
        @test findfirst(CompatibleSearchQuery(dna""), seq) === 1:0
        @test findfirst(CompatibleSearchQuery(dna"AC"), seq) === 1:2
        @test findnext(CompatibleSearchQuery(dna"AC"), seq, 2) === 5:6
        #@test findfirst(CompatibleSearchQuery(dna"AC"), seq, 2, 5) === nothing
        @test findfirst(CompatibleSearchQuery(dna"TG"), seq) === nothing
        @test findfirst(CompatibleSearchQuery(dna"TN"), seq) === 4:5
        @test findfirst(CompatibleSearchQuery(dna"ACG"), seq)  === 1:3
        @test findnext(CompatibleSearchQuery(dna"ACG"), seq, 2) === 5:7
        @test findfirst(CompatibleSearchQuery(seq), seq) === 1:lastindex(seq)

        @test findfirst(CompatibleSearchQuery(dna""), dna"") === 1:0
        @test findnext(CompatibleSearchQuery(dna""), dna"", -1) === 1:0
        @test findnext(CompatibleSearchQuery(dna""), dna"", 2) === nothing

        @test first(findfirst(CompatibleSearchQuery(dna""), seq)) === 1
        @test first(findfirst(CompatibleSearchQuery(dna"AC"), seq)) === 1
        @test first(findnext(CompatibleSearchQuery(dna"AC"), seq, 2)) === 5
        #@test findnext(CompatibleSearchQuery(dna"AC"), seq, 2) === nothing

        query = CompatibleSearchQuery(dna"ACG")
        @test findfirst(query, seq) === 1:3
        @test findnext(query, seq, 2) === 5:7
        #@test findfirst(query, seq, 2, 6) === nothing
        @test first(findfirst(query, seq)) === 1
        @test first(findnext(query, seq, 2)) === 5
        #@test findfirst(query, seq, 2, 6) === nothing
    end

    @testset "backward" begin
        @test findlast(CompatibleSearchQuery(dna""), seq) === 8:7
        @test findlast(CompatibleSearchQuery(dna"AC"), seq) === 5:6
        @test findprev(CompatibleSearchQuery(dna"AC"), seq, 5) === 1:2
        #@test findlast(CompatibleSearchQuery(dna"AC"), seq, 5, 2) === nothing
        @test findlast(CompatibleSearchQuery(dna"TG"), seq) === nothing
        @test findlast(CompatibleSearchQuery(dna"TN"), seq) === 4:5
        @test findlast(CompatibleSearchQuery(dna"ACG"), seq) === 5:7
        @test findprev(CompatibleSearchQuery(dna"ACG"), seq, 6) === 1:3
        @test findlast(CompatibleSearchQuery(seq), seq) === 1:lastindex(seq)

        @test findlast(CompatibleSearchQuery(dna""), dna"") === 1:0
        @test findprev(CompatibleSearchQuery(dna""), dna"", 2) === 1:0
        @test findprev(CompatibleSearchQuery(dna""), dna"", -1) === nothing

        @test first(findlast(CompatibleSearchQuery(dna""), seq)) === 8
        @test first(findlast(CompatibleSearchQuery(dna"AC"), seq)) === 5
        @test first(findprev(CompatibleSearchQuery(dna"AC"), seq, 5)) === 1
        #@test findlast(CompatibleSearchQuery(dna"AC"), seq, 5, 2) === nothing

        query = CompatibleSearchQuery(dna"ACG")
        @test findlast(query, seq) === 5:7
        @test findprev(query, seq, 6) === 1:3
        #@test findlast(query, seq, 2, 6) === nothing
        @test first(findlast(query, seq)) === 5
        @test first(findprev(query, seq, 6)) === 1
        #@test findlast(query, seq, 6, 2) === nothing
    end

    @testset "occursin" begin
        @test occursin(CompatibleSearchQuery(dna"ACG"), dna"GGGTACACGTTT") == true
        @test occursin(CompatibleSearchQuery(dna"TGT"), dna"GGGTACACGTGT") == true
        @test occursin(CompatibleSearchQuery(dna"GGG"), dna"GGGTACACGTGT") == true
        @test occursin(CompatibleSearchQuery(dna"AAA"), dna"GGGTACACGTGT") == false
    end

end
