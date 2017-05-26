@testset "Order" begin
    @test DNAKmer("AA") < DNAKmer("AC") < DNAKmer("AG") < DNAKmer("AT") < DNAKmer("CA")
    @test RNAKmer("AA") < RNAKmer("AC") < RNAKmer("AG") < RNAKmer("AU") < RNAKmer("CA")
end
