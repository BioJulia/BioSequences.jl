@testset "De Bruijn Neighbors" begin
    @test collect(neighbors(DNAKmer("ACG")))  ==
        DNAKmer["CGA",  "CGC",  "CGG",  "CGT" ]
    @test collect(neighbors(DNAKmer("GGGG"))) ==
        DNAKmer["GGGA", "GGGC", "GGGG", "GGGT"]
    @test collect(neighbors(RNAKmer("ACG")))  ==
        RNAKmer["CGA",  "CGC",  "CGG",  "CGU" ]
    @test collect(neighbors(RNAKmer("GGGG"))) ==
        RNAKmer["GGGA", "GGGC", "GGGG", "GGGU"]
end
