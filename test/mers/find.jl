@testset "Find" begin
    kmer = DNAMer("ACGAG")
    bigkmer = BigDNAMer("ACGAG")

    @test findnext(isequal(DNA_A), kmer, 1) == 1
    @test findnext(isequal(DNA_C), kmer, 1) == 2
    @test findnext(isequal(DNA_G), kmer, 1) == 3
    @test findnext(isequal(DNA_T), kmer, 1) === nothing
    @test findnext(isequal(DNA_A), kmer, 2) == 4
    
    @test findnext(isequal(DNA_A), bigkmer, 1) == 1
    @test findnext(isequal(DNA_C), bigkmer, 1) == 2
    @test findnext(isequal(DNA_G), bigkmer, 1) == 3
    @test findnext(isequal(DNA_T), bigkmer, 1) === nothing
    @test findnext(isequal(DNA_A), bigkmer, 2) == 4

    @test_throws BoundsError findnext(isequal(DNA_A), kmer, 0)
    @test findnext(isequal(DNA_A), kmer, 6) === nothing
    
    @test_throws BoundsError findnext(isequal(DNA_A), bigkmer, 0)
    @test findnext(isequal(DNA_A), bigkmer, 6) === nothing

    @test findprev(isequal(DNA_A), kmer, 5) == 4
    @test findprev(isequal(DNA_C), kmer, 5) == 2
    @test findprev(isequal(DNA_G), kmer, 5) == 5
    @test findprev(isequal(DNA_T), kmer, 5) === nothing
    @test findprev(isequal(DNA_G), kmer, 4) == 3
    
    @test findprev(isequal(DNA_A), bigkmer, 5) == 4
    @test findprev(isequal(DNA_C), bigkmer, 5) == 2
    @test findprev(isequal(DNA_G), bigkmer, 5) == 5
    @test findprev(isequal(DNA_T), bigkmer, 5) === nothing
    @test findprev(isequal(DNA_G), bigkmer, 4) == 3

    @test findprev(isequal(DNA_A), kmer, 0) === nothing
    @test_throws BoundsError findprev(isequal(DNA_A), kmer, 6)
    
    @test findprev(isequal(DNA_A), bigkmer, 0) === nothing
    @test_throws BoundsError findprev(isequal(DNA_A), bigkmer, 6)

    @test findfirst(isequal(DNA_A), kmer) == 1
    @test findfirst(isequal(DNA_G), kmer) == 3
    @test findlast(isequal(DNA_A), kmer) == 4
    @test findlast(isequal(DNA_G), kmer) == 5
    
    @test findfirst(isequal(DNA_A), bigkmer) == 1
    @test findfirst(isequal(DNA_G), bigkmer) == 3
    @test findlast(isequal(DNA_A), bigkmer) == 4
    @test findlast(isequal(DNA_G), bigkmer) == 5
end
