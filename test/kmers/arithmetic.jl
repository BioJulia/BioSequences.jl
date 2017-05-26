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
end
