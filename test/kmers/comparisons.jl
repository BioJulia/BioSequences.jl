@testset "Comparisons" begin
    @testset "Equality" begin
        function check_seq_kmer_equality(len)
            a = DNAMer(random_dna_kmer(len))
            b = DNASequence(a)
            return a == b && b == a
        end

        for len in [1, 10, 32]
            @test all(Bool[check_seq_kmer_equality(len) for _ in 1:reps])
        end

        # True negatives
        @test DNAMer("ACG") != RNAMer("ACG")
        @test DNAMer("T")   != RNAMer("U")
        @test DNAMer("AC")  != DNAMer("AG")
        @test RNAMer("AC")  != RNAMer("AG")

        @test DNAMer("ACG") != rna"ACG"
        @test DNAMer("T")   != rna"U"
        @test DNAMer("AC")  != dna"AG"
        @test RNAMer("AC")  != rna"AG"

        @test rna"ACG" != DNAMer("ACG")
        @test rna"U"   != DNAMer("T")
        @test dna"AG"  != DNAMer("AC")
        @test rna"AG"  != RNAMer("AC")
    end

    @testset "Inequality" begin
        for len in [1, 10, 32]
            @test  isless(DNAMer{1}(UInt64(0)), DNAMer{1}(UInt64(1)))
            @test !isless(DNAMer{1}(UInt64(0)), DNAMer{1}(UInt64(0)))
            @test !isless(DNAMer{1}(UInt64(1)), DNAMer{1}(UInt64(0)))

            @test  isless(RNAMer{1}(UInt64(0)), RNAMer{1}(UInt64(1)))
            @test !isless(RNAMer{1}(UInt64(0)), RNAMer{1}(UInt64(0)))
            @test !isless(RNAMer{1}(UInt64(1)), RNAMer{1}(UInt64(0)))
        end
    end

    @testset "Hash" begin
        kmers = map(DNAMer, ["AAAA", "AACT", "ACGT", "TGCA"])
        for x in kmers, y in kmers
            @test (x == y) == (hash(x) == hash(y))
        end
        kmers = map(RNAMer, ["AAAA", "AACU", "ACGU", "UGCA"])
        for x in kmers, y in kmers
            @test (x == y) == (hash(x) == hash(y))
        end
    end
end
