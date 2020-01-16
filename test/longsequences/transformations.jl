@testset "Transformations" begin
    function test_reverse(A, seq)
        revseq = reverse(LongSequence{A}(seq))
        @test convert(String, revseq) == reverse(seq)
    end

    function test_dna_complement(A, seq)
        comp = complement(LongSequence{A}(seq))
        @test convert(String, comp) == dna_complement(seq)
    end

    function test_rna_complement(A, seq)
        comp = complement(LongSequence{A}(seq))
        @test convert(String, comp) == rna_complement(seq)
    end

    function test_dna_revcomp(A, seq)
        revcomp = reverse_complement(LongSequence{A}(seq))
        @test convert(String, revcomp) == reverse(dna_complement(seq))
    end

    function test_rna_revcomp(A, seq)
        revcomp = reverse_complement(LongSequence{A}(seq))
        @test convert(String, revcomp) == reverse(rna_complement(seq))
    end

    @testset "Reverse" begin
        for len in [0, 1, 10, 32, 1000, 10000, 100000], _ in 1:10
            test_reverse(DNAAlphabet{4}, random_dna(len))
            test_reverse(RNAAlphabet{4}, random_rna(len))
            test_reverse(AminoAcidAlphabet, random_aa(len))
            test_reverse(CharAlphabet, random_aa(len))

            probs = [0.25, 0.25, 0.25, 0.25, 0.00]
            test_reverse(DNAAlphabet{2}, random_dna(len, probs))
            test_reverse(RNAAlphabet{2}, random_rna(len, probs))
        end
    end

    @testset "Complement" begin
        for len in [0, 1, 10, 32, 1000, 10000, 100000], _ in 1:10
            test_dna_complement(DNAAlphabet{4}, random_dna(len))
            test_rna_complement(RNAAlphabet{4}, random_rna(len))

            probs = [0.25, 0.25, 0.25, 0.25, 0.00]
            test_dna_complement(DNAAlphabet{2}, random_dna(len, probs))
            test_rna_complement(RNAAlphabet{2}, random_rna(len, probs))
        end
        seq_string = join(rand("-ACGTSWKMYRBDHVN", 1000))
        seq = complement(LongSequence{DNAAlphabet{4}}(seq_string))
        @test String(seq) == dna_complement(seq_string)

        seq_string = join(rand("-ACGUSWKMYRBDHVN", 1000))
        seq = complement(LongSequence{RNAAlphabet{4}}(seq_string))
        @test String(seq) == rna_complement(seq_string)
    end

    @testset "Reverse complement" begin
        for len in [0, 1, 10, 32, 1000, 10000, 100000], _ in 1:10
            test_dna_revcomp(DNAAlphabet{4}, random_dna(len))
            test_rna_revcomp(RNAAlphabet{4}, random_rna(len))

            probs = [0.25, 0.25, 0.25, 0.25, 0.00]
            test_dna_revcomp(DNAAlphabet{2}, random_dna(len, probs))
            test_rna_revcomp(RNAAlphabet{2}, random_rna(len, probs))
        end
        seq_string = join(rand("-ACGTSWKMYRBDHVN", 1000))
        seq = reverse_complement(LongSequence{DNAAlphabet{4}}(seq_string))
        @test String(seq) == reverse(dna_complement(seq_string))

        seq_string = join(rand("-ACGUSWKMYRBDHVN", 1000))
        seq = reverse_complement(LongSequence{RNAAlphabet{4}}(seq_string))
        @test String(seq) == reverse(rna_complement(seq_string))
    end

    @testset "Map" begin
        seq = dna""
        @test map(identity, seq) == dna""
        seq = dna"AAA"
        @test map(x -> DNA_C, seq) == dna"CCC"
        seq = dna"ACGTNACGTN"
        @test map(x -> x == DNA_N ? DNA_A : x, seq) == dna"ACGTAACGTA"
        @test seq == dna"ACGTNACGTN"
        @test map!(x -> x == DNA_N ? DNA_A : x, seq) === seq
        @test seq == dna"ACGTAACGTA"

        seq = rna""
        @test map(identity, seq) == rna""
        seq = rna"AAA"
        @test map(x -> RNA_C, seq) == rna"CCC"
        seq = rna"ACGUNACGUN"
        @test map(x -> x == RNA_N ? RNA_A : x, seq) == rna"ACGUAACGUA"
        @test seq == rna"ACGUNACGUN"
        @test map!(x -> x == RNA_N ? RNA_A : x, seq) === seq
        @test seq == rna"ACGUAACGUA"

        seq = aa""
        @test map(identity, seq) == aa""
        seq = aa"MMM"
        @test map(x -> AA_P, seq) == aa"PPP"
        seq = aa"XRNDCQXE"
        @test map(x -> x == AA_X ? AA_A : x, seq) == aa"ARNDCQAE"
        @test seq == aa"XRNDCQXE"
        @test map!(x -> x == AA_X ? AA_A : x, seq) === seq
        @test seq == aa"ARNDCQAE"
    end

    @testset "Filter" begin
        seq = dna""
        @test filter(x -> true, seq) == dna""
        @test filter(x -> false, seq) == dna""
        seq = dna"AAA"
        @test filter(x -> x == DNA_A, seq) == dna"AAA"
        @test filter(x -> x == DNA_C, seq) == dna""
        seq = dna"ACGTNACGTN"
        @test filter(x -> x == DNA_N, seq) == dna"NN"
        @test filter(x -> x != DNA_N, seq) == dna"ACGTACGT"
        @test seq == dna"ACGTNACGTN"
        @test filter!(x -> x != DNA_N, seq) == seq
        @test seq == dna"ACGTACGT"

        for len in 1:50, _ in 1:10
            str = random_dna(len)
            seq = LongDNASeq(str)
            @test filter(x -> x == DNA_N, seq) ==
                LongDNASeq(filter(x -> x == 'N', str))
            @test filter(x -> x != DNA_N, seq) ==
                LongDNASeq(filter(x -> x != 'N', str))
        end

        seq = rna""
        @test filter(x -> true, seq) == rna""
        @test filter(x -> false, seq) == rna""
        seq = rna"AAA"
        @test filter(x -> x == RNA_A, seq) == rna"AAA"
        @test filter(x -> x == RNA_C, seq) == rna""
        seq = rna"ACGUNACGUN"
        @test filter(x -> x == RNA_N, seq) == rna"NN"
        @test filter(x -> x != RNA_N, seq) == rna"ACGUACGU"
        @test seq == rna"ACGUNACGUN"
        @test filter!(x -> x != RNA_N, seq) == seq
        @test seq == rna"ACGUACGU"

        seq = aa""
        @test filter(x -> true, seq) == aa""
        @test filter(x -> false, seq) == aa""
        seq = aa"PPP"
        @test filter(x -> x == AA_P, seq) == aa"PPP"
        @test filter(x -> x != AA_P, seq) == aa""
        seq = aa"ARNDCQXGHILKMFPXTWYVOUX"
        @test filter(x -> x == AA_X, seq) == aa"XXX"
        @test filter(x -> x != AA_X, seq) == aa"ARNDCQGHILKMFPTWYVOU"
        @test seq == aa"ARNDCQXGHILKMFPXTWYVOUX"
        @test filter!(x -> x != AA_X, seq) == seq
        @test seq == aa"ARNDCQGHILKMFPTWYVOU"
    end
end
