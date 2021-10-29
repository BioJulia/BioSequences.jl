@testset "Approximate" begin
    seq = dna"ACGTACG"

    @testset "forward" begin
        ApproximateSearchQuery(dna"", 0)
        @test findnext(ApproximateSearchQuery(dna"", 0), seq) === 1:0
        @test findnext(ApproximateSearchQuery(dna"AC", 0), seq) === 1:2
        @test findnext(ApproximateSearchQuery(dna"AC", 0), seq, 2) === 5:6
        #@test approxsearch(seq, dna"AC", 0, 2, 5) === 0:-1 make as view
        @test findnext(ApproximateSearchQuery(dna"AT", 0), seq) === 0:-1
        @test findnext(ApproximateSearchQuery(dna"AT", 1), seq) === 1:1
        @test findnext(ApproximateSearchQuery(dna"AT", 1), seq, 2) === 3:4
        #@test approxsearch(seq, dna"AT", 1, 2, 3) === 0:-1 make as view
        @test findnext(ApproximateSearchQuery(dna"NG", 0), seq) === 2:3
        @test findnext(ApproximateSearchQuery(dna"NG", 1), seq) === 1:1
        @test findnext(ApproximateSearchQuery(dna"GN", 0), seq) === 3:4
        @test findnext(ApproximateSearchQuery(dna"GN", 1), seq) === 1:1
        
        @test findnext(ApproximateSearchQuery(dna"ACG", 0), seq) === 1:3
        @test findnext(ApproximateSearchQuery(dna"ACG", 1), seq) === 1:2
        @test findnext(ApproximateSearchQuery(dna"ACG", 2), seq) === 1:1
        @test findnext(ApproximateSearchQuery(dna"ACG", 3), seq) === 1:0
        @test findnext(ApproximateSearchQuery(dna"ACG", 4), seq) === 1:0
    end

    @testset "backward" begin
        # TODO: maybe this should return 8:7 like rsearch
        @test approxrsearch(seq, dna"", 0) === 7:6
        @test approxrsearch(seq, dna"AC", 0) === 5:6
        @test approxrsearch(seq, dna"AC", 0, 5) === 1:2
        @test approxrsearch(seq, dna"AC", 0, 5, 2) === nothing
        @test approxrsearch(seq, dna"AT", 0) === nothing
        @test approxrsearch(seq, dna"AT", 1) === 5:6
        @test approxrsearch(seq, dna"AT", 1, 6) === 5:6
        @test approxrsearch(seq, dna"AT", 1, 5) === 5:5
        @test approxrsearch(seq, dna"AT", 1, 3, 2) === nothing
        @test approxrsearch(seq, dna"NG", 0) === 6:7
        @test approxrsearch(seq, dna"NG", 0, 6) === 2:3
        @test approxrsearch(seq, dna"GN", 0) === 3:4
        @test approxrsearch(seq, dna"GN", 1) === 7:7
        @test approxrsearch(seq, dna"ACG", 0) === 5:7
        @test approxrsearch(seq, dna"ACG", 1) === 6:7
        @test approxrsearch(seq, dna"ACG", 2) === 7:7
        @test approxrsearch(seq, dna"ACG", 3) === nothing
        @test approxrsearch(seq, dna"ACG", 4) === nothing

        # TODO: maybe this should return 8 like rsearchindex
        @test approxrsearchindex(seq, dna"", 0) === 7
        @test approxrsearchindex(seq, dna"AC", 0) === 5
        @test approxrsearchindex(seq, dna"AC", 0, 5) === 1
        @test approxrsearchindex(seq, dna"AC", 0, 5, 2) === nothing

        query = ApproximateSearchQuery(dna"ACG")
        @test approxrsearch(seq, query, 1, 7) === 6:7
        @test approxrsearch(seq, query, 1, 6) === 5:6
        @test approxrsearch(seq, query, 1, 6, 6) === nothing
        @test approxrsearchindex(seq, query, 1, 7) === 6
        @test approxrsearchindex(seq, query, 1, 6) === 5
        @test approxrsearchindex(seq, query, 1, 6, 6) === nothing
    end
end
