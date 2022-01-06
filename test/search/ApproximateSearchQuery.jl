@testset "Approximate" begin
    seq = dna"ACGTACG"

    @testset "iscompatible" begin
        @testset "forward" begin
            @test findnext(ApproximateSearchQuery(dna"", iscompatible), 0, seq, 1) === 1:0
            @test findnext(ApproximateSearchQuery(dna"AC", iscompatible), 0, seq, 1) === 1:2
            @test findnext(ApproximateSearchQuery(dna"AC", iscompatible), 0, seq, 2) === 5:6
            
            @test findnext(ApproximateSearchQuery(dna"AT", iscompatible), 0, seq, 1) === nothing
            @test findnext(ApproximateSearchQuery(dna"AT", iscompatible), 1, seq, 1) === 1:1
            @test findnext(ApproximateSearchQuery(dna"AT", iscompatible), 1, seq, 2) === 3:4
            
            @test findnext(ApproximateSearchQuery(dna"NG", iscompatible), 0, seq, 1) === 2:3
            @test findnext(ApproximateSearchQuery(dna"NG", iscompatible), 1, seq, 1) === 1:1
            @test findnext(ApproximateSearchQuery(dna"GN", iscompatible), 0, seq, 1) === 3:4
            @test findnext(ApproximateSearchQuery(dna"GN", iscompatible), 1, seq, 1) === 1:1
            
            @test findnext(ApproximateSearchQuery(dna"ACG", iscompatible), 0, seq, 1) === 1:3
            @test findnext(ApproximateSearchQuery(dna"ACG", iscompatible), 1, seq, 1) === 1:2
            @test findnext(ApproximateSearchQuery(dna"ACG", iscompatible), 2, seq, 1) === 1:1
            @test findnext(ApproximateSearchQuery(dna"ACG", iscompatible), 3, seq, 1) === 1:0
            @test findnext(ApproximateSearchQuery(dna"ACG", iscompatible), 4, seq, 1) === 1:0
        end
        
        @testset "backward" begin
            # TODO: maybe this should return 8:7 like rsearch
            
            @test findprev(ApproximateSearchQuery(dna"", iscompatible), 0, seq, 7) === 7:6
            @test findprev(ApproximateSearchQuery(dna"AC", iscompatible), 0, seq, 7) === 5:6
            @test findprev(ApproximateSearchQuery(dna"AC", iscompatible), 0, seq, 5) === 1:2
            
            @test findprev(ApproximateSearchQuery(dna"AT", iscompatible), 0, seq, 7) === nothing
            @test findprev(ApproximateSearchQuery(dna"AT", iscompatible), 1, seq, 7) === 5:6
            @test findprev(ApproximateSearchQuery(dna"AT", iscompatible), 1, seq, 6) === 5:6
            @test findprev(ApproximateSearchQuery(dna"AT", iscompatible), 1, seq, 5) === 5:5
            
            @test findprev(ApproximateSearchQuery(dna"NG", iscompatible), 0, seq, 7) === 6:7
            @test findprev(ApproximateSearchQuery(dna"NG", iscompatible), 0, seq, 6) === 2:3
            @test findprev(ApproximateSearchQuery(dna"GN", iscompatible), 0, seq, 7) === 3:4
            @test findprev(ApproximateSearchQuery(dna"GN", iscompatible), 1, seq, 7) === 7:7
            
            @test findprev(ApproximateSearchQuery(dna"ACG", iscompatible), 0, seq, 7) === 5:7
            @test findprev(ApproximateSearchQuery(dna"ACG", iscompatible), 1, seq, 7) === 6:7
            @test findprev(ApproximateSearchQuery(dna"ACG", iscompatible), 2, seq, 7) === 7:7
            @test findprev(ApproximateSearchQuery(dna"ACG", iscompatible), 3, seq, 7) === 7:6
            @test findprev(ApproximateSearchQuery(dna"ACG", iscompatible), 4, seq, 7) === 7:6
        end
    end

    @testset "isequal" begin
        @testset "forward" begin
            @test findnext(ApproximateSearchQuery(dna""), 0, seq, 1) === 1:0
            @test findnext(ApproximateSearchQuery(dna"AC"), 0, seq, 1) === 1:2
            @test findnext(ApproximateSearchQuery(dna"AC"), 0, seq, 2) === 5:6
            
            @test findnext(ApproximateSearchQuery(dna"AT"), 0, seq, 1) === nothing
            @test findnext(ApproximateSearchQuery(dna"AT"), 1, seq, 1) === 1:1
            @test findnext(ApproximateSearchQuery(dna"AT"), 1, seq, 2) === 3:4
            
            @test findnext(ApproximateSearchQuery(dna"NG"), 0, seq, 1) === nothing
            @test findnext(ApproximateSearchQuery(dna"NG"), 1, seq, 1) === 2:3
            @test findnext(ApproximateSearchQuery(dna"GN"), 0, seq, 1) === nothing
            @test findnext(ApproximateSearchQuery(dna"GN"), 1, seq, 1) === 3:3
            
            @test findnext(ApproximateSearchQuery(dna"ACG"), 0, seq, 1) === 1:3
            @test findnext(ApproximateSearchQuery(dna"ACG"), 1, seq, 1) === 1:2
            @test findnext(ApproximateSearchQuery(dna"ACG"), 2, seq, 1) === 1:1
            @test findnext(ApproximateSearchQuery(dna"ACG"), 3, seq, 1) === 1:0
            @test findnext(ApproximateSearchQuery(dna"ACG"), 4, seq, 1) === 1:0
        end
        
        @testset "backward" begin
            # TODO: maybe this should return 8:7 like rsearch
            @test findprev(ApproximateSearchQuery(dna""), 0, seq, 7) === 7:6
            @test findprev(ApproximateSearchQuery(dna"AC"), 0, seq, 7) === 5:6
            @test findprev(ApproximateSearchQuery(dna"AC"), 0, seq, 5) === 1:2
            
            @test findprev(ApproximateSearchQuery(dna"AT"), 0, seq, 7) === nothing
            @test findprev(ApproximateSearchQuery(dna"AT"), 1, seq, 7) === 5:6
            @test findprev(ApproximateSearchQuery(dna"AT"), 1, seq, 6) === 5:6
            @test findprev(ApproximateSearchQuery(dna"AT"), 1, seq, 5) === 5:5
            
            @test findprev(ApproximateSearchQuery(dna"NG"), 0, seq, 7) === nothing
            @test findprev(ApproximateSearchQuery(dna"NG"), 1, seq, 6) === 3:3
            @test findprev(ApproximateSearchQuery(dna"GN"), 0, seq, 7) === nothing
            @test findprev(ApproximateSearchQuery(dna"GN"), 1, seq, 7) === 7:7
            
            @test findprev(ApproximateSearchQuery(dna"ACG"), 0, seq, 7) === 5:7
            @test findprev(ApproximateSearchQuery(dna"ACG"), 1, seq, 7) === 6:7
            @test findprev(ApproximateSearchQuery(dna"ACG"), 2, seq, 7) === 7:7
            @test findprev(ApproximateSearchQuery(dna"ACG"), 3, seq, 7) === 7:6
            @test findprev(ApproximateSearchQuery(dna"ACG"), 4, seq, 7) === 7:6
        end
    end
end
