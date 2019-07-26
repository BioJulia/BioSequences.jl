@testset "De Bruijn Neighbors" begin
    @test collect(neighbors(DNAMer("ACG")))  == map(DNAMer, ["CGA",  "CGC",  "CGG",  "CGT"])
    @test collect(neighbors(DNAMer("GGGG"))) == map(DNAMer, ["GGGA", "GGGC", "GGGG", "GGGT"])
    @test collect(neighbors(RNAMer("ACG")))  == map(RNAMer, ["CGA",  "CGC",  "CGG",  "CGU"])
    @test collect(neighbors(RNAMer("GGGG"))) == map(RNAMer, ["GGGA", "GGGC", "GGGG", "GGGU"])
    
    @test collect(neighbors(BigDNAMer("ACG")))  == map(BigDNAMer, ["CGA",  "CGC",  "CGG",  "CGT"])
    @test collect(neighbors(BigDNAMer("GGGG"))) == map(BigDNAMer, ["GGGA", "GGGC", "GGGG", "GGGT"])
    @test collect(neighbors(BigRNAMer("ACG")))  == map(BigRNAMer, ["CGA",  "CGC",  "CGG",  "CGU"])
    @test collect(neighbors(BigRNAMer("GGGG"))) == map(BigRNAMer, ["GGGA", "GGGC", "GGGG", "GGGU"])
end
