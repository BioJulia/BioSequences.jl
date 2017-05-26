@testset "Exact" begin
    seq = dna"ACGTACG"

    @testset "forward" begin
        @test search(seq, DNA_A) === 1:1
        @test search(seq, DNA_N) === 1:1
        @test search(seq, DNA_T) === 4:4
        @test search(seq, DNA_T, 5) === 0:-1
        @test search(seq, DNA_A, 2, 6) === 5:5

        @test search(seq, dna"") === 1:0
        @test search(seq, dna"AC") === 1:2
        @test search(seq, dna"AC", 2) === 5:6
        @test search(seq, dna"AC", 2, 5) === 0:-1
        @test search(seq, dna"TG") === 0:-1
        @test search(seq, dna"TN") === 4:5
        @test search(seq, dna"ACG") === 1:3
        @test search(seq, dna"ACG", 2) === 5:7
        @test search(seq, seq) === 1:endof(seq)

        @test search(dna"", dna"") === 1:0
        @test search(dna"", dna"", -1) === 1:0
        @test search(dna"", dna"", 2) === 0:-1

        @test searchindex(seq, dna"") === 1
        @test searchindex(seq, dna"AC") === 1
        @test searchindex(seq, dna"AC", 2) === 5
        @test searchindex(seq, dna"AC", 2, 5) === 0

        query = ExactSearchQuery(dna"ACG")
        @test search(seq, query) === 1:3
        @test search(seq, query, 2) === 5:7
        @test search(seq, query, 2, 6) === 0:-1
        @test searchindex(seq, query) === 1
        @test searchindex(seq, query, 2) === 5
        @test searchindex(seq, query, 2, 6) === 0
    end

    @testset "backward" begin
        @test rsearch(seq, DNA_A) === 5:5
        @test rsearch(seq, DNA_N) === 7:7
        @test rsearch(seq, DNA_T) === 4:4
        @test rsearch(seq, DNA_T, 3) === 0:-1
        @test rsearch(seq, DNA_T, 6, 3) === 4:4

        @test rsearch(seq, dna"") === 8:7
        @test rsearch(seq, dna"AC") === 5:6
        @test rsearch(seq, dna"AC", 5) === 1:2
        @test rsearch(seq, dna"AC", 5, 2) === 0:-1
        @test rsearch(seq, dna"TG") === 0:-1
        @test rsearch(seq, dna"TN") === 4:5
        @test rsearch(seq, dna"ACG") === 5:7
        @test rsearch(seq, dna"ACG", 6) === 1:3
        @test rsearch(seq, seq) === 1:endof(seq)

        @test rsearch(dna"", dna"") === 1:0
        @test rsearch(dna"", dna"", 2) === 1:0
        @test rsearch(dna"", dna"", -1) === 0:-1

        @test rsearchindex(seq, dna"") === 8
        @test rsearchindex(seq, dna"AC") === 5
        @test rsearchindex(seq, dna"AC", 5) === 1
        @test rsearchindex(seq, dna"AC", 5, 2) === 0

        query = ExactSearchQuery(dna"ACG")
        @test rsearch(seq, query) === 5:7
        @test rsearch(seq, query, 6) === 1:3
        @test rsearch(seq, query, 2, 6) === 0:-1
        @test rsearchindex(seq, query) === 5
        @test rsearchindex(seq, query, 6) === 1
        @test rsearchindex(seq, query, 6, 2) === 0
    end
end
