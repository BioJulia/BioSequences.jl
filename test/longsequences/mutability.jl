@testset "Mutability" begin
    @testset "setindex!" begin
        s = dna"ACGT"
        s[1] = DNA_A
        @test s == dna"ACGT"
        s[1] = DNA_C
        @test s == dna"CCGT"
        s[2] = DNA_T
        @test s == dna"CTGT"
        s[4] = DNA_A
        @test s == dna"CTGA"
        @test_throws BoundsError s[0]
        @test_throws BoundsError s[5]

        s = dna"ACGTACGT"
        s[3:5] = DNA_A
        @test s == dna"ACAAACGT"
        s[3:5] = DNA_T
        @test s == dna"ACTTTCGT"
        s[1:8] = DNA_C
        @test s == dna"CCCCCCCC"
        s[:] = DNA_G
        @test s == dna"GGGGGGGG"
        @test_throws BoundsError s[0:3]
        @test_throws BoundsError s[5:10]

        s = dna"ACGTACGT"
        s[[1,2,3]] = DNA_T
        @test s == dna"TTTTACGT"
        s[[1,5,8]] = DNA_C
        @test s == dna"CTTTCCGC"
        @test_throws BoundsError s[[3,9]] = DNA_A

        s = dna"ACGT"
        s[[true, false, false, true]] = DNA_G
        @test s == dna"GCGG"
        s[trues(4)] = DNA_A
        @test s == dna"AAAA"
        @test_throws BoundsError s[[true, false, false]] = DNA_G
        @test_throws BoundsError s[[true, false, false, false, true]] = DNA_G

        s = dna"ACGTACGT"
        s[2:3] = dna"AA"
        @test s == dna"AAATACGT"
        s[7:8] = dna"CC"
        @test s == dna"AAATACCC"
        s[:] = dna"AACCGGTT"
        @test s == dna"AACCGGTT"
        @test_throws BoundsError       s[0:1] = dna"AA"
        @test_throws DimensionMismatch s[3:4] = dna"A"

        s = dna"ACGTACGT"
        s[[1,4]] = dna"TA"
        @test s == dna"TCGAACGT"
        s[[2,3,5]] = dna"CAT"
        @test s == dna"TCAATCGT"
        @test_throws BoundsError       s[[1,2,9]] = dna"AAA"
        @test_throws DimensionMismatch s[[1,2,8]] = dna"AA"

        s = dna"ACGT"
        s[[true,false,true,false]] = dna"TT"
        @test s == dna"TCTT"
        s[trues(4)] = dna"AAAA"
        @test s == dna"AAAA"
        @test_throws BoundsError       s[[true,false,true]] = dna"TT"
        @test_throws DimensionMismatch s[[true,false,true,true]] = dna"TT"
    end

    @testset "resize!" begin
        seq = dna""
        resize!(seq, 100)
        @test length(seq) == 100
        resize!(seq, 200)
        @test length(seq) == 200
        resize!(seq,  10)
        @test length(seq) == 10
        @test_throws ArgumentError resize!(seq, -1)
    end

    @testset "empty!" begin
        seq = dna"ACG"
        @test empty!(seq) == dna""
        @test length(seq) == 0
    end

    @testset "push!" begin
        seq = dna""
        @test push!(seq, DNA_A) == dna"A"
        @test push!(seq, DNA_C) == dna"AC"
        @test seq == dna"AC"
    end

    @testset "pushfirst!" begin
        seq = dna""
        @test pushfirst!(seq, DNA_A) == dna"A"
        @test pushfirst!(seq, DNA_C) == dna"CA"
        @test seq == dna"CA"
    end

    @testset "pop!" begin
        seq = dna"ACGT"
        @test pop!(seq) === DNA_T
        @test seq == dna"ACG"
        @test pop!(seq) === DNA_G
        @test seq == dna"AC"
        @test_throws ArgumentError pop!(dna"")
    end

    @testset "popfirst!" begin
        seq = dna"ACGT"
        @test popfirst!(seq) === DNA_A
        @test seq == dna"CGT"
        @test popfirst!(seq) === DNA_C
        @test seq == dna"GT"
        @test_throws ArgumentError popfirst!(dna"")
    end

    @testset "insert!" begin
        seq = dna"ACGT"
        @test insert!(seq, 2, DNA_G) == dna"AGCGT"
        @test insert!(seq, 5, DNA_A) == dna"AGCGAT"
        @test_throws BoundsError insert!(seq, 10, DNA_T)
    end

    @testset "deleteat!" begin
        seq = dna"ACGT"
        @test deleteat!(seq, 1) == dna"CGT"
        @test deleteat!(seq, 2) == dna"CT"
        @test_throws BoundsError deleteat!(seq, 10)

        seq = dna"ACGTACGT"
        @test deleteat!(seq, 3:5) == dna"ACCGT"
        @test_throws BoundsError deleteat!(seq, 10:12)
    end

    @testset "append!" begin
        seq = dna""
        @test append!(seq, dna"A") == dna"A"
        @test append!(seq, dna"ACG") == dna"AACG"
    end

    @testset "copyto!" begin
        seq = dna"GGG"
        @test copyto!(seq, dna"ACG") == dna"ACG"
        @test copyto!(seq, dna"TTA") == dna"TTA"

        seq = dna"TCCC"
        @test copyto!(seq, 2, dna"TT", 1) == dna"TTTC"
        seq = dna"TCCC"
        @test copyto!(seq, 2, dna"TT", 1, 1) == dna"TTCC"

        seq = dna"ACGT"
        @test copyto!(seq, seq) == dna"ACGT"
        @test copyto!(seq, 1, seq, 3, 2) == dna"GTGT"
        seq = dna"ACGT"
        @test copyto!(seq, 3, seq, 1, 2) == dna"ACAC"
    end

    @testset "orphan!" begin
        seq = repeat(dna"ACGT", 8)
        subseq = seq[16:17]
        BioSequences.orphan!(subseq)
        @test subseq == dna"TA"
    end
end
