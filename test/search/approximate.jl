@testset "Approximate" begin
    seq = dna"ACGTACG"

    @testset "forward" begin
        @test approxsearch(seq, dna"", 0) === 1:0
        @test approxsearch(seq, dna"AC", 0) === 1:2
        @test approxsearch(seq, dna"AC", 0, 2) === 5:6
        @test approxsearch(seq, dna"AC", 0, 2, 5) === 0:-1
        @test approxsearch(seq, dna"AT", 0) === 0:-1
        @test approxsearch(seq, dna"AT", 1) === 1:1
        @test approxsearch(seq, dna"AT", 1, 2) === 3:4
        @test approxsearch(seq, dna"AT", 1, 2, 3) === 0:-1
        @test approxsearch(seq, dna"NG", 0) === 2:3
        @test approxsearch(seq, dna"NG", 1) === 1:1
        @test approxsearch(seq, dna"GN", 0) === 3:4
        @test approxsearch(seq, dna"GN", 1) === 1:1
        @test approxsearch(seq, dna"ACG", 0) === 1:3
        @test approxsearch(seq, dna"ACG", 1) === 1:2
        @test approxsearch(seq, dna"ACG", 2) === 1:1
        @test approxsearch(seq, dna"ACG", 3) === 1:0
        @test approxsearch(seq, dna"ACG", 4) === 1:0

        @test approxsearchindex(seq, dna"", 0) === 1
        @test approxsearchindex(seq, dna"AC", 0) === 1
        @test approxsearchindex(seq, dna"AC", 0, 2) === 5
        @test approxsearchindex(seq, dna"AC", 0, 2, 5) === 0

        query = ApproximateSearchQuery(dna"ACG")
        @test approxsearch(seq, query, 1) === 1:2
        @test approxsearch(seq, query, 1, 2) === 2:3
        @test approxsearch(seq, query, 1, 2, 2) === 0:-1
        @test approxsearchindex(seq, query, 1) === 1
        @test approxsearchindex(seq, query, 1, 2) === 2
        @test approxsearchindex(seq, query, 1, 2, 2) === 0
    end

    @testset "backward" begin
        # TODO: maybe this should return 8:7 like rsearch
        @test approxrsearch(seq, dna"", 0) === 7:6
        @test approxrsearch(seq, dna"AC", 0) === 5:6
        @test approxrsearch(seq, dna"AC", 0, 5) === 1:2
        @test approxrsearch(seq, dna"AC", 0, 5, 2) === 0:-1
        @test approxrsearch(seq, dna"AT", 0) === 0:-1
        @test approxrsearch(seq, dna"AT", 1) === 5:6
        @test approxrsearch(seq, dna"AT", 1, 6) === 5:6
        @test approxrsearch(seq, dna"AT", 1, 5) === 5:5
        @test approxrsearch(seq, dna"AT", 1, 3, 2) === 0:-1
        @test approxrsearch(seq, dna"NG", 0) === 6:7
        @test approxrsearch(seq, dna"NG", 0, 6) === 2:3
        @test approxrsearch(seq, dna"GN", 0) === 3:4
        @test approxrsearch(seq, dna"GN", 1) === 7:7
        @test approxrsearch(seq, dna"ACG", 0) === 5:7
        @test approxrsearch(seq, dna"ACG", 1) === 6:7
        @test approxrsearch(seq, dna"ACG", 2) === 7:7
        @test approxrsearch(seq, dna"ACG", 3) === 7:6
        @test approxrsearch(seq, dna"ACG", 4) === 7:6

        # TODO: maybe this should return 8 like rsearchindex
        @test approxrsearchindex(seq, dna"", 0) === 7
        @test approxrsearchindex(seq, dna"AC", 0) === 5
        @test approxrsearchindex(seq, dna"AC", 0, 5) === 1
        @test approxrsearchindex(seq, dna"AC", 0, 5, 2) === 0

        query = ApproximateSearchQuery(dna"ACG")
        @test approxrsearch(seq, query, 1, 7) === 6:7
        @test approxrsearch(seq, query, 1, 6) === 5:6
        @test approxrsearch(seq, query, 1, 6, 6) === 0:-1
        @test approxrsearchindex(seq, query, 1, 7) === 6
        @test approxrsearchindex(seq, query, 1, 6) === 5
        @test approxrsearchindex(seq, query, 1, 6, 6) === 0
    end
end
