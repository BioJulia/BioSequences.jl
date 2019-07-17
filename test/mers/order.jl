@testset "Order" begin
    @test DNAMer("AA") < DNAMer("AC") < DNAMer("AG") < DNAMer("AT") < DNAMer("CA")
    @test RNAMer("AA") < RNAMer("AC") < RNAMer("AG") < RNAMer("AU") < RNAMer("CA")
end
