@testset "De Bruijn Neighbors" begin
    @test collect(neighbors(DNAKmer("ACG")))  ==
        map(DNAKmer, ["CGA",  "CGC",  "CGG",  "CGT"])
    @test collect(neighbors(DNAKmer("GGGG"))) ==
        map(DNAKmer, ["GGGA", "GGGC", "GGGG", "GGGT"])
    @test collect(neighbors(RNAKmer("ACG")))  ==
        map(RNAKmer, ["CGA",  "CGC",  "CGG",  "CGU"])
    @test collect(neighbors(RNAKmer("GGGG"))) ==
        map(RNAKmer, ["GGGA", "GGGC", "GGGG", "GGGU"])
end
