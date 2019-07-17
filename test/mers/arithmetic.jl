@testset "Arithmetic" begin
    x = DNAMer("AA")
    @test x - 1 == x + (-1) == x - 0x01 == DNAMer("TT")
    @test x + 1 == x - (-1) == x + 0x01 == DNAMer("AC")

    x = DNAMer("TT")
    @test x - 1 == x + (-1) == x - 0x01 == DNAMer("TG")
    @test x + 1 == x - (-1) == x + 0x01 == DNAMer("AA")

    base = DNAMer("AAA")
    offset = 0
    nucs = "ACGT"
    for a in nucs, b in nucs, c in nucs
        @test base + offset == DNAMer(string(a, b, c))
        offset += 1
    end

    @testset "typemin" begin
        for k in 1:32
            @test typemin(DNAMer{k}) === DNAMer("A"^k)
            @test typemin(RNAMer{k}) === RNAMer("A"^k)
        end
        @test_throws ArgumentError typemin(DNAMer{0})
        @test_throws ArgumentError typemin(DNAMer{33})
    end

    @testset "typemax" begin
        for k in 1:32
            @test typemax(DNAMer{k}) === DNAMer("T"^k)
            @test typemax(RNAMer{k}) === RNAMer("U"^k)
        end
        @test_throws ArgumentError typemax(DNAMer{0})
        @test_throws ArgumentError typemax(DNAMer{33})
    end
end
