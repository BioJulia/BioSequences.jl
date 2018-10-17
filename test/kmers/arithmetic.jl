@testset "Arithmetic" begin
    x = DNAKmer("AA")
    @test x - 1 == x + (-1) == x - 0x01 == DNAKmer("TT")
    @test x + 1 == x - (-1) == x + 0x01 == DNAKmer("AC")

    x = DNAKmer("TT")
    @test x - 1 == x + (-1) == x - 0x01 == DNAKmer("TG")
    @test x + 1 == x - (-1) == x + 0x01 == DNAKmer("AA")

    base = DNAKmer("AAA")
    offset = 0
    nucs = "ACGT"
    for a in nucs, b in nucs, c in nucs
        @test base + offset == DNAKmer(string(a, b, c))
        offset += 1
    end

    @testset "typemin" begin
        for k in 1:32
            @test typemin(DNAKmer{k}) === DNAKmer("A"^k)
            @test typemin(RNAKmer{k}) === RNAKmer("A"^k)
        end
        @test_throws ArgumentError typemin(DNAKmer{0})
        @test_throws ArgumentError typemin(DNAKmer{33})
    end

    @testset "typemax" begin
        for k in 1:32
            @test typemax(DNAKmer{k}) === DNAKmer("T"^k)
            @test typemax(RNAKmer{k}) === RNAKmer("U"^k)
        end
        @test_throws ArgumentError typemax(DNAKmer{0})
        @test_throws ArgumentError typemax(DNAKmer{33})
    end
end
