@testset "Subsequence construction" begin
    function test_subseq(A, seq)
        bioseq = GeneralSequence{A}(seq)
        for _ in 1:100
            part = random_interval(1, lastindex(seq))
            @test convert(String, bioseq[part]) == seq[part]
        end
    end

    for len in [1, 10, 32, 1000, 10000, 100000]
        test_subseq(DNAAlphabet{4}, random_dna(len))
        test_subseq(RNAAlphabet{4}, random_rna(len))
        test_subseq(AminoAcidAlphabet, random_aa(len))

        probs = [0.25, 0.25, 0.25, 0.25, 0.00]
        test_subseq(DNAAlphabet{2}, random_dna(len, probs))
        test_subseq(RNAAlphabet{2}, random_rna(len, probs))
    end

    # Subsequence from range
    @test RNASequence(rna"AUCGAUCG", 5:8) == RNASequence("AUCG")
    @test DNASequence(dna"ATCGATCG", 5:8) == DNASequence("ATCG")

    # Invalid ranges
    @test_throws Exception RNASequence(rna"AUCGAUCG", 5:10)
    @test_throws Exception DNASequence(dna"ATCGATCG", 5:10)

    # Empty ranges
    @test RNASequence(rna"AUCGAUCG", 5:4) == RNASequence()
    @test DNASequence(dna"ATCGATCG", 5:4) == DNASequence()

    # Subsequence of subsequence
    @test dna"ACGTAG"[4:end][1:2] == dna"TA"
    @test dna"ACGTAG"[4:end][2:3] == dna"AG"
    @test_throws Exception dna"ACGTAG"[4:end][0:1]
    @test_throws Exception dna"ACGTAG"[4:end][3:4]
end
