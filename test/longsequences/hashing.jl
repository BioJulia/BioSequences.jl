@testset "Hash" begin
    seq = dna"ACGTACGT"
    @test isa(hash(seq), UInt64)
    @test hash(seq) === hash(dna"ACGTACGT")
    @test hash(seq) !== hash(seq[1:6])
    @test hash(seq) !== hash(seq[1:7])
    @test hash(seq) === hash(seq[1:8])
    @test hash(seq[1:4]) === hash(dna"ACGT")
    @test hash(seq[2:4]) === hash(dna"CGT")
    @test hash(seq[3:4]) === hash(dna"GT")
    @test hash(seq[4:4]) === hash(dna"T")
    @test hash(seq[5:8]) === hash(dna"ACGT")

    @test hash(dna"") !== hash(dna"A")
    @test hash(dna"A") !== hash(dna"AA")
    @test hash(dna"AA") !== hash(dna"AAA")
    @test hash(dna"AAA") !== hash(dna"AAAA")

    for n in 1:200, seq in [dna"A", dna"AC", dna"ACG", dna"ACGT", dna"ACGTN"]
        @test hash(seq^n) === hash((dna""     * seq^n)[1:end])
        @test hash(seq^n) === hash((dna"T"    * seq^n)[2:end])
        @test hash(seq^n) === hash((dna"TT"   * seq^n)[3:end])
        @test hash(seq^n) === hash((dna"TTT"  * seq^n)[4:end])
        @test hash(seq^n) === hash((dna"TTTT" * seq^n)[5:end])

        @test hash(seq^n) === hash((seq^n * dna""    )[1:end  ])
        @test hash(seq^n) === hash((seq^n * dna"T"   )[1:end-1])
        @test hash(seq^n) === hash((seq^n * dna"TT"  )[1:end-2])
        @test hash(seq^n) === hash((seq^n * dna"TTT" )[1:end-3])
        @test hash(seq^n) === hash((seq^n * dna"TTTT")[1:end-4])
    end

    @test hash(rna"AAUU") === hash(rna"AAUU")
    @test hash(rna"AAUUAA"[3:5]) === hash(rna"UUA")
    @test hash(aa"MTTQAPMFTQPLQ") === hash(aa"MTTQAPMFTQPLQ")
    @test hash(aa"MTTQAPMFTQPLQ"[5:10]) === hash(aa"APMFTQ")

    @testset "MinHash" begin
        seq = DNASequence(random_dna(1000))
        h = minhash(seq, 10, 100)

        @test length(h) == 100
        @test h == minhash(seq, 10, 100)

        @test_throws BoundsError h[101]
    end

end
