@testset "GC content" begin
    @test gc_content(dna"") === 0.0
    @test gc_content(dna"AATA") === 0.0
    @test gc_content(dna"ACGT") === 0.5
    @test gc_content(dna"CGGC") === 1.0
    @test gc_content(dna"ACATTGTGTATAACAAAAGG") === 6 / 20
    @test gc_content(dna"GAGGCGTTTATCATC"[2:end]) === 6 / 14

    @test gc_content(DNAKmer("AATA")) === 0.0
    @test gc_content(DNAKmer("ACGT")) === 0.5
    @test gc_content(DNAKmer("CGGC")) === 1.0
    @test gc_content(DNAKmer("ACATTGTGTATAACAAAAGG")) === 6 / 20

    @test gc_content(rna"") === 0.0
    @test gc_content(rna"AAUA") === 0.0
    @test gc_content(rna"ACGU") === 0.5
    @test gc_content(rna"CGGC") === 1.0
    @test gc_content(rna"ACAUUGUGUAUAACAAAAGG") === 6 / 20
    @test gc_content(rna"GAGGCGUUUAUCAUC"[2:end]) === 6 / 14

    @test gc_content(RNAKmer("AAUA")) === 0.0
    @test gc_content(RNAKmer("ACGU")) === 0.5
    @test gc_content(RNAKmer("CGGC")) === 1.0
    @test gc_content(RNAKmer("ACAUUGUGUAUAACAAAAGG")) === 6 / 20

    @test_throws Exception gc_content(aa"ARN")

    Random.seed!(1234)
    for _ in 1:200
        s = randdnaseq(rand(1:100))
        @test gc_content(s) === count(isGC, s) / length(s)
        @test gc_content(GeneralSequence{DNAAlphabet{2}}(s)) === count(isGC, s) / length(s)
        i = rand(1:endof(s))
        j = rand(i-1:endof(s))
        @test gc_content(s[i:j]) === (j < i ? 0.0 : count(isGC, s[i:j]) / (j - i + 1))
    end
end
