@testset "De Bruijn Neighbors" begin
    @test collect(neighbors(DNAMer("ACG")))  ==
        map(DNAMer, ["CGA",  "CGC",  "CGG",  "CGT"])
    @test collect(neighbors(DNAMer("GGGG"))) ==
        map(DNAMer, ["GGGA", "GGGC", "GGGG", "GGGT"])
    @test collect(neighbors(RNAMer("ACG")))  ==
        map(RNAMer, ["CGA",  "CGC",  "CGG",  "CGU"])
    @test collect(neighbors(RNAMer("GGGG"))) ==
        map(RNAMer, ["GGGA", "GGGC", "GGGG", "GGGU"])
end
