@testset "Approximate" begin
    seq = dna"ACGTACG"

    @testset "forward" begin
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
        
        @test findprev(ApproximateSearchQuery(dna"", 0), seq) === 7:6
        @test findprev(ApproximateSearchQuery(dna"AC", 0), seq) === 5:6
        @test findprev(ApproximateSearchQuery(dna"AC", 0), seq, 5) === 1:2
        #@test approxsearch(seq, dna"AC", 0, 5, 2) === 0:-1 make as view
        @test findprev(ApproximateSearchQuery(dna"AT", 0), seq) === 0:-1
        @test findprev(ApproximateSearchQuery(dna"AT", 1), seq) === 5:6
        @test findprev(ApproximateSearchQuery(dna"AT", 1), seq, 6) === 5:6
        @test findprev(ApproximateSearchQuery(dna"AT", 1), seq, 5) === 5:5
        #@test approxsearch(seq, dna"AT", 1, 3, 2) === 0:-1 make as view
        @test findprev(ApproximateSearchQuery(dna"NG", 0), seq) === 6:7
        @test findprev(ApproximateSearchQuery(dna"NG", 0), seq, 6) === 2:3
        @test findprev(ApproximateSearchQuery(dna"GN", 0), seq) === 3:4
        @test findprev(ApproximateSearchQuery(dna"GN", 1), seq) === 7:7
        
        @test findprev(ApproximateSearchQuery(dna"ACG", 0), seq) === 5:7
        @test findprev(ApproximateSearchQuery(dna"ACG", 1), seq) === 6:7
        @test findprev(ApproximateSearchQuery(dna"ACG", 2), seq) === 7:7
        @test findprev(ApproximateSearchQuery(dna"ACG", 3), seq) === 7:6
        @test findprev(ApproximateSearchQuery(dna"ACG", 4), seq) === 7:6
    end
end
