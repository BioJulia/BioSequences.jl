@testset "Subsequence construction" begin
    function test_subseq(A, seq)
        bioseq = LongSequence{A}(seq)
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
    @test LongRNASeq(rna"AUCGAUCG", 5:8) == LongRNASeq("AUCG")
    @test LongDNASeq(dna"ATCGATCG", 5:8) == LongDNASeq("ATCG")

    # Invalid ranges
    @test_throws Exception LongRNASeq(rna"AUCGAUCG", 5:10)
    @test_throws Exception LongDNASeq(dna"ATCGATCG", 5:10)

    # Empty ranges
    @test LongRNASeq(rna"AUCGAUCG", 5:4) == LongRNASeq()
    @test LongDNASeq(dna"ATCGATCG", 5:4) == LongDNASeq()

    # Subsequence of subsequence
    @test dna"ACGTAG"[4:end][1:2] == dna"TA"
    @test dna"ACGTAG"[4:end][2:3] == dna"AG"
    @test_throws Exception dna"ACGTAG"[4:end][0:1]
    @test_throws Exception dna"ACGTAG"[4:end][3:4]
end
