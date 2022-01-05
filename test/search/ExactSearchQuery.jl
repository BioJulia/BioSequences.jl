@testset "ExactSearchQuery" begin
    seq = dna"ACGTACG"
    
    @testset "isequal" begin
        @testset "forward" begin
            @test findfirst(ExactSearchQuery(dna""), seq) === 1:0
            @test findfirst(ExactSearchQuery(dna"AC"), seq) === 1:2
            @test findnext(ExactSearchQuery(dna"AC"), seq, 2) === 5:6
            #@test findfirst(ExactSearchQuery(dna"AC"), seq, 2, 5) === nothing
            @test findfirst(ExactSearchQuery(dna"TG"), seq) === nothing
            @test findfirst(ExactSearchQuery(dna"TN"), seq) === nothing
            @test findfirst(ExactSearchQuery(dna"ACG"), seq)  === 1:3
            @test findnext(ExactSearchQuery(dna"ACG"), seq, 2) === 5:7
            @test findfirst(ExactSearchQuery(seq), seq) === 1:lastindex(seq)

            @test findfirst(ExactSearchQuery(dna""), dna"") === 1:0
            @test findnext(ExactSearchQuery(dna""), dna"", -1) === 1:0
            @test findnext(ExactSearchQuery(dna""), dna"", 2) === nothing

            @test first(findfirst(ExactSearchQuery(dna""), seq)) === 1
            @test first(findfirst(ExactSearchQuery(dna"AC"), seq)) === 1
            @test first(findnext(ExactSearchQuery(dna"AC"), seq, 2)) === 5
            #@test findnext(ExactSearchQuery(dna"AC"), seq, 2) === nothing

            query = ExactSearchQuery(dna"ACG")
            @test findfirst(query, seq) === 1:3
            @test findnext(query, seq, 2) === 5:7
            #@test findfirst(query, seq, 2, 6) === nothing
            @test first(findfirst(query, seq)) === 1
            @test first(findnext(query, seq, 2)) === 5
            #@test findfirst(query, seq, 2, 6) === nothing
        end

        @testset "backward" begin
            @test findlast(ExactSearchQuery(dna""), seq) === 8:7
            @test findlast(ExactSearchQuery(dna"AC"), seq) === 5:6
            @test findprev(ExactSearchQuery(dna"AC"), seq, 5) === 1:2
            #@test findlast(ExactSearchQuery(dna"AC"), seq, 5, 2) === nothing
            @test findlast(ExactSearchQuery(dna"TG"), seq) === nothing
            @test findlast(ExactSearchQuery(dna"TN"), seq) === nothing
            @test findlast(ExactSearchQuery(dna"ACG"), seq) === 5:7
            @test findprev(ExactSearchQuery(dna"ACG"), seq, 6) === 1:3
            @test findlast(ExactSearchQuery(seq), seq) === 1:lastindex(seq)

            @test findlast(ExactSearchQuery(dna""), dna"") === 1:0
            @test findprev(ExactSearchQuery(dna""), dna"", 2) === 1:0
            @test findprev(ExactSearchQuery(dna""), dna"", -1) === nothing

            @test first(findlast(ExactSearchQuery(dna""), seq)) === 8
            @test first(findlast(ExactSearchQuery(dna"AC"), seq)) === 5
            @test first(findprev(ExactSearchQuery(dna"AC"), seq, 5)) === 1
            #@test findlast(ExactSearchQuery(dna"AC"), seq, 5, 2) === nothing

            query = ExactSearchQuery(dna"ACG")
            @test findlast(query, seq) === 5:7
            @test findprev(query, seq, 6) === 1:3
            #@test findlast(query, seq, 2, 6) === nothing
            @test first(findlast(query, seq)) === 5
            @test first(findprev(query, seq, 6)) === 1
            #@test findlast(query, seq, 6, 2) === nothing
        end

        @testset "occursin" begin
            @test occursin(ExactSearchQuery(dna"ACG"), dna"GGGTACACGTTT") == true
            @test occursin(ExactSearchQuery(dna"TGT"), dna"GGGTACACGTGT") == true
            @test occursin(ExactSearchQuery(dna"GGG"), dna"GGGTACACGTGT") == true
            @test occursin(ExactSearchQuery(dna"AAA"), dna"GGGTACACGTGT") == false
        end
    end
    
    @testset "iscompatible" begin
        
        @testset "forward" begin
            @test findfirst(ExactSearchQuery(dna"", iscompatible), seq) === 1:0
            @test findfirst(ExactSearchQuery(dna"AC", iscompatible), seq) === 1:2
            @test findnext(ExactSearchQuery(dna"AC", iscompatible), seq, 2) === 5:6
            #@test findfirst(ExactSearchQuery(dna"AC"), seq, 2, 5) === nothing
            @test findfirst(ExactSearchQuery(dna"TG", iscompatible), seq) === nothing
            @test findfirst(ExactSearchQuery(dna"TN", iscompatible), seq) === 4:5
            @test findfirst(ExactSearchQuery(dna"ACG", iscompatible), seq)  === 1:3
            @test findnext(ExactSearchQuery(dna"ACG", iscompatible), seq, 2) === 5:7
            @test findfirst(ExactSearchQuery(seq, iscompatible), seq) === 1:lastindex(seq)

            @test findfirst(ExactSearchQuery(dna"", iscompatible), dna"") === 1:0
            @test findnext(ExactSearchQuery(dna"", iscompatible), dna"", -1) === 1:0
            @test findnext(ExactSearchQuery(dna"", iscompatible), dna"", 2) === nothing

            @test first(findfirst(ExactSearchQuery(dna"", iscompatible), seq)) === 1
            @test first(findfirst(ExactSearchQuery(dna"AC", iscompatible), seq)) === 1
            @test first(findnext(ExactSearchQuery(dna"AC", iscompatible), seq, 2)) === 5
            #@test findnext(ExactSearchQuery(dna"AC"), seq, 2) === nothing

            query = ExactSearchQuery(dna"ACG", iscompatible)
            @test findfirst(query, seq) === 1:3
            @test findnext(query, seq, 2) === 5:7
            #@test findfirst(query, seq, 2, 6) === nothing
            @test first(findfirst(query, seq)) === 1
            @test first(findnext(query, seq, 2)) === 5
            #@test findfirst(query, seq, 2, 6) === nothing
        end

        @testset "backward" begin
            @test findlast(ExactSearchQuery(dna"", iscompatible), seq) === 8:7
            @test findlast(ExactSearchQuery(dna"AC", iscompatible), seq) === 5:6
            @test findprev(ExactSearchQuery(dna"AC", iscompatible), seq, 5) === 1:2
            #@test findlast(ExactSearchQuery(dna"AC"), seq, 5, 2) === nothing
            @test findlast(ExactSearchQuery(dna"TG", iscompatible), seq) === nothing
            @test findlast(ExactSearchQuery(dna"TN", iscompatible), seq) === 4:5
            @test findlast(ExactSearchQuery(dna"ACG", iscompatible), seq) === 5:7
            @test findprev(ExactSearchQuery(dna"ACG", iscompatible), seq, 6) === 1:3
            @test findlast(ExactSearchQuery(seq, iscompatible), seq) === 1:lastindex(seq)

            @test findlast(ExactSearchQuery(dna"", iscompatible), dna"") === 1:0
            @test findprev(ExactSearchQuery(dna"", iscompatible), dna"", 2) === 1:0
            @test findprev(ExactSearchQuery(dna"", iscompatible), dna"", -1) === nothing

            @test first(findlast(ExactSearchQuery(dna"", iscompatible), seq)) === 8
            @test first(findlast(ExactSearchQuery(dna"AC", iscompatible), seq)) === 5
            @test first(findprev(ExactSearchQuery(dna"AC", iscompatible), seq, 5)) === 1
            #@test findlast(ExactSearchQuery(dna"AC"), seq, 5, 2) === nothing

            query = ExactSearchQuery(dna"ACG", iscompatible)
            @test findlast(query, seq) === 5:7
            @test findprev(query, seq, 6) === 1:3
            #@test findlast(query, seq, 2, 6) === nothing
            @test first(findlast(query, seq)) === 5
            @test first(findprev(query, seq, 6)) === 1
            #@test findlast(query, seq, 6, 2) === nothing
        end

        @testset "occursin" begin
            @test occursin(ExactSearchQuery(dna"ACG", iscompatible), dna"GGGTACACGTTT") == true
            @test occursin(ExactSearchQuery(dna"TGT", iscompatible), dna"GGGTACACGTGT") == true
            @test occursin(ExactSearchQuery(dna"GGG", iscompatible), dna"GGGTACACGTGT") == true
            @test occursin(ExactSearchQuery(dna"AAA", iscompatible), dna"GGGTACACGTGT") == false
        end
    end

end
