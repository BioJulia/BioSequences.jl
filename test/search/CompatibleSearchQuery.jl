@testset "CompatibleQuery" begin
    seq = dna"ACGTACG"
    
    @testset "forward" begin
        @test findfirst(CompatibleQuery(dna""), seq) === 1:0
        @test findfirst(CompatibleQuery(dna"AC"), seq) === 1:2
        @test findnext(CompatibleQuery(dna"AC"), seq, 2) === 5:6
        #@test findfirst(CompatibleQuery(dna"AC"), seq, 2, 5) === nothing
        @test findfirst(CompatibleQuery(dna"TG"), seq) === nothing
        @test findfirst(CompatibleQuery(dna"TN"), seq) === 4:5
        @test findfirst(CompatibleQuery(dna"ACG"), seq)  === 1:3
        @test findnext(CompatibleQuery(dna"ACG"), seq, 2) === 5:7
        @test findfirst(CompatibleQuery(seq), seq) === 1:lastindex(seq)

        @test findfirst(CompatibleQuery(dna""), dna"") === 1:0
        @test findnext(CompatibleQuery(dna""), dna"", -1) === 1:0
        @test findnext(CompatibleQuery(dna""), dna"", 2) === nothing

        @test first(findfirst(CompatibleQuery(dna""), seq)) === 1
        @test first(findfirst(CompatibleQuery(dna"AC"), seq)) === 1
        @test first(findnext(CompatibleQuery(dna"AC"), seq, 2)) === 5
        #@test findnext(CompatibleQuery(dna"AC"), seq, 2) === nothing

        query = CompatibleQuery(dna"ACG")
        @test findfirst(query, seq) === 1:3
        @test findnext(query, seq, 2) === 5:7
        #@test findfirst(query, seq, 2, 6) === nothing
        @test first(findfirst(query, seq)) === 1
        @test first(findnext(query, seq, 2)) === 5
        #@test findfirst(query, seq, 2, 6) === nothing
    end

    @testset "backward" begin
        @test findlast(CompatibleQuery(dna""), seq) === 8:7
        @test findlast(CompatibleQuery(dna"AC"), seq) === 5:6
        @test findprev(CompatibleQuery(dna"AC"), seq, 5) === 1:2
        #@test findlast(CompatibleQuery(dna"AC"), seq, 5, 2) === nothing
        @test findlast(CompatibleQuery(dna"TG"), seq) === nothing
        @test findlast(CompatibleQuery(dna"TN"), seq) === 4:5
        @test findlast(CompatibleQuery(dna"ACG"), seq) === 5:7
        @test findprev(CompatibleQuery(dna"ACG"), seq, 6) === 1:3
        @test findlast(CompatibleQuery(seq), seq) === 1:lastindex(seq)

        @test findlast(CompatibleQuery(dna""), dna"") === 1:0
        @test findprev(CompatibleQuery(dna""), dna"", 2) === 1:0
        @test findprev(CompatibleQuery(dna""), dna"", -1) === nothing

        @test first(findlast(CompatibleQuery(dna""), seq)) === 8
        @test first(findlast(CompatibleQuery(dna"AC"), seq)) === 5
        @test first(findprev(CompatibleQuery(dna"AC"), seq, 5)) === 1
        #@test findlast(CompatibleQuery(dna"AC"), seq, 5, 2) === nothing

        query = CompatibleQuery(dna"ACG")
        @test findlast(query, seq) === 5:7
        @test findprev(query, seq, 6) === 1:3
        #@test findlast(query, seq, 2, 6) === nothing
        @test first(findlast(query, seq)) === 5
        @test first(findprev(query, seq, 6)) === 1
        #@test findlast(query, seq, 6, 2) === nothing
    end

    @testset "occursin" begin
        @test occursin(CompatibleQuery(dna"ACG"), dna"GGGTACACGTTT") == true
        @test occursin(CompatibleQuery(dna"TGT"), dna"GGGTACACGTGT") == true
        @test occursin(CompatibleQuery(dna"GGG"), dna"GGGTACACGTGT") == true
        @test occursin(CompatibleQuery(dna"AAA"), dna"GGGTACACGTGT") == false
    end

end
