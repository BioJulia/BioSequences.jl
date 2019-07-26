@testset "Order" begin
    @test DNAMer("AA") < DNAMer("AC") < DNAMer("AG") < DNAMer("AT") < DNAMer("CA")
    @test RNAMer("AA") < RNAMer("AC") < RNAMer("AG") < RNAMer("AU") < RNAMer("CA")
    
    @test BigDNAMer("AA") < BigDNAMer("AC") < BigDNAMer("AG") < BigDNAMer("AT") < BigDNAMer("CA")
    @test BigRNAMer("AA") < BigRNAMer("AC") < BigRNAMer("AG") < BigRNAMer("AU") < BigRNAMer("CA")
end
