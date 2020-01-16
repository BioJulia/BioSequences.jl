@testset "Copy" begin
    function test_copy(A, seq)
        @test convert(String, copy(LongSequence{A}(seq))) == seq
    end

    for len in [1, 10, 16, 32, 1000, 10000]
        test_copy(DNAAlphabet{4}, random_dna(len))
        test_copy(RNAAlphabet{4}, random_rna(len))
        test_copy(AminoAcidAlphabet, random_aa(len))

        probs = [0.25, 0.25, 0.25, 0.25, 0.00]
        test_copy(DNAAlphabet{2}, random_dna(len, probs))
        test_copy(RNAAlphabet{2}, random_rna(len, probs))
    end

    seq = dna"ACGTACGTACGT"
    subseq = seq[3:6]
    @test copy(subseq) == dna"GTAC"
end

@testset "Copy!" begin
    function test_copy!(seq, src)
        @test String(src) == String(copy!(copy(seq), src))
    end
    # Needed because conversion to String truncates vector.
    function test_copy!(seq, src::Vector)
        @test String(copy(src)) == String(copy!(copy(seq), src))
    end

    probs = [0.25, 0.25, 0.25, 0.25, 0.00]
    dna2 = LongSequence{DNAAlphabet{2}}(6)
    dna4 = LongSequence{DNAAlphabet{4}}(6)
    rna2 = LongSequence{RNAAlphabet{2}}(6)
    rna4 = LongSequence{RNAAlphabet{4}}(6)
    aa = LongSequence{AminoAcidAlphabet}(6)
    charseq = LongSequence{CharAlphabet}(6)
    for dtype in [Vector{UInt8}, Vector{Char}, String, Test.GenericString]
        for len in [0, 1, 10, 16, 32, 100, 5]
            test_copy!(dna2, dtype(random_dna(len, probs)))
            test_copy!(dna4, dtype(random_dna(len)))
            test_copy!(rna2, dtype(random_rna(len, probs)))
            test_copy!(rna4, dtype(random_rna(len)))
            test_copy!(aa, dtype(random_aa(len)))
            test_copy!(charseq, dtype(random_aa(len)))
        end
    end
    test_copy!(charseq, "ϐʌ⨝W")
end

@testset "Concatenation" begin
    function test_concatenation(A, chunks)
        parts = UnitRange{Int}[]
        for i in 1:lastindex(chunks)
            start = rand(1:length(chunks[i]))
            stop = rand(start:length(chunks[i]))
            push!(parts, start:stop)
        end
        str = string([chunk[parts[i]] for (i, chunk) in enumerate(chunks)]...)
        seq = *([LongSequence{A}(chunk)[parts[i]] for (i, chunk) in enumerate(chunks)]...)
        @test convert(String, seq) == uppercase(str)
    end

    for _ in 1:100
        n = rand(1:10)
        chunks = [random_dna(rand(1:100)) for _ in 1:n]
        test_concatenation(DNAAlphabet{4}, chunks)

        chunks = [random_rna(rand(1:100)) for _ in 1:n]
        test_concatenation(RNAAlphabet{4}, chunks)

        chunks = [random_aa(rand(1:100)) for _ in 1:n]
        test_concatenation(AminoAcidAlphabet, chunks)

        probs = [0.25, 0.25, 0.25, 0.25, 0.00]
        chunks = [random_dna(rand(1:100), probs) for _ in 1:n]
        test_concatenation(DNAAlphabet{2}, chunks)

        chunks = [random_rna(rand(1:100), probs) for _ in 1:n]
        test_concatenation(RNAAlphabet{2}, chunks)
    end
end

@testset "Repetition" begin
    function test_repetition(A, chunk)
        start = rand(1:length(chunk))
        stop = rand(start:length(chunk))
        n = rand(1:10)
        str = chunk[start:stop] ^ n
        seq = LongSequence{A}(chunk)[start:stop] ^ n
        @test convert(String, seq) == uppercase(str)
    end

    for _ in 1:10
        chunk = random_dna(rand(1:100))
        test_repetition(DNAAlphabet{4}, chunk)

        chunk = random_rna(rand(1:100))
        test_repetition(RNAAlphabet{4}, chunk)

        chunk = random_aa(rand(1:100))
        test_repetition(AminoAcidAlphabet, chunk)

        probs = [0.25, 0.25, 0.25, 0.25, 0.00]
        chunk = random_dna(rand(1:100), probs)
        test_repetition(DNAAlphabet{2}, chunk)

        chunk = random_rna(rand(1:100), probs)
        test_repetition(RNAAlphabet{2}, chunk)
    end
end

@testset "Length" begin
    for len in [0, 1, 2, 3, 10, 16, 32, 1000, 10000]
        seq = LongDNASeq(random_dna(len))
        @test length(seq) === lastindex(seq) === len

        seq = LongRNASeq(random_rna(len))
        @test length(seq) === lastindex(seq) === len

        seq = LongAminoAcidSeq(random_aa(len))
        @test length(seq) === lastindex(seq) === len
    end

    @test length(char"いろはabc") === 6
end

@testset "Access" begin
    dna_seq = dna"ACTG"

    @test dna_seq[1] === DNA_A
    @test dna_seq[2] === DNA_C
    @test dna_seq[3] === DNA_T
    @test dna_seq[4] === DNA_G

    # Access indexes out of bounds
    @test_throws BoundsError dna_seq[-1]
    @test_throws BoundsError dna_seq[0]
    @test_throws BoundsError dna_seq[5]

    @test dna"ACTGNACTGN"[1:5] == dna"ACTGN"
    @test dna"ACTGNACTGN"[5:1] == dna""

    rna_seq = rna"ACUG"
    @test rna_seq[1] === RNA_A
    @test rna_seq[2] === RNA_C
    @test rna_seq[3] === RNA_U
    @test rna_seq[4] === RNA_G

    # Access indexes out of bounds
    @test_throws BoundsError rna_seq[-1]
    @test_throws BoundsError rna_seq[0]
    @test_throws BoundsError rna_seq[5]

    @test rna"ACUGNACUGN"[1:5] == rna"ACUGN"
    @test rna"ACUGNACUGN"[5:1] == rna""

    @test aa"KSAAV"[3] == AA_A
    @test char"いろはにほ"[3] == 'は'
end

@testset "Equality" begin
    @testset "DNA" begin
        a = dna"ACTGN"
        b = dna"ACTGN"
        @test a == b
        @test dna"ACTGN" == dna"ACTGN"
        @test dna"ACTGN" != dna"ACTGA"
        @test dna"ACTGN" != dna"ACTG"
        @test dna"ACTG"  != dna"ACTGN"

        a = dna"ACGTNACGTN"
        b = dna"""
        ACGTN
        ACGTN
        """
        @test a == b
    end

    @testset "RNA" begin
        a = rna"ACUGN"
        b = rna"ACUGN"
        @test a == b
        @test rna"ACUGN" == rna"ACUGN"
        @test rna"ACUGN" != rna"ACUGA"
        @test rna"ACUGN" != rna"ACUG"
        @test rna"ACUG"  != rna"ACUGN"

        a = rna"ACUGNACUGN"
        b = rna"""
        ACUGN
        ACUGN
        """
        @test a == b
    end

    @testset "AminoAcid" begin
        a = aa"ARNDCQEGHILKMFPSTWYVX"
        b = aa"ARNDCQEGHILKMFPSTWYVX"
        @test a == b
        @test a != aa"ARNDCQEGHILKMFPSTWYXV"
        @test a != aa"ARNDCQEGHLKMFPSTWYVX"

        b = aa"""
        ARNDCQEGHI
        LKMFPSTWYV
        X
        """
        @test a == b
    end
end
