@testset "Comparisons" begin
    @testset "Equality" begin
        function check_seq_kmer_equality(len)
            a = DNAKmer(random_dna_kmer(len))
            b = DNASequence(a)
            return a == b && b == a
        end

        for len in [1, 10, 32]
            @test all(Bool[check_seq_kmer_equality(len) for _ in 1:reps])
        end

        # True negatives
        @test DNAKmer("ACG") != RNAKmer("ACG")
        @test DNAKmer("T")   != RNAKmer("U")
        @test DNAKmer("AC")  != DNAKmer("AG")
        @test RNAKmer("AC")  != RNAKmer("AG")

        @test DNAKmer("ACG") != rna"ACG"
        @test DNAKmer("T")   != rna"U"
        @test DNAKmer("AC")  != dna"AG"
        @test RNAKmer("AC")  != rna"AG"

        @test rna"ACG" != DNAKmer("ACG")
        @test rna"U"   != DNAKmer("T")
        @test dna"AG"  != DNAKmer("AC")
        @test rna"AG"  != RNAKmer("AC")
    end

    @testset "Inequality" begin
        for len in [1, 10, 32]
            @test  isless(DNAKmer{1}(UInt64(0)), DNAKmer{1}(UInt64(1)))
            @test !isless(DNAKmer{1}(UInt64(0)), DNAKmer{1}(UInt64(0)))
            @test !isless(DNAKmer{1}(UInt64(1)), DNAKmer{1}(UInt64(0)))

            @test  isless(RNAKmer{1}(UInt64(0)), RNAKmer{1}(UInt64(1)))
            @test !isless(RNAKmer{1}(UInt64(0)), RNAKmer{1}(UInt64(0)))
            @test !isless(RNAKmer{1}(UInt64(1)), RNAKmer{1}(UInt64(0)))
        end
    end

    @testset "Hash" begin
        kmers = map(DNAKmer, ["AAAA", "AACT", "ACGT", "TGCA"])
        for x in kmers, y in kmers
            @test (x == y) == (hash(x) == hash(y))
        end
        kmers = map(RNAKmer, ["AAAA", "AACU", "ACGU", "UGCA"])
        for x in kmers, y in kmers
            @test (x == y) == (hash(x) == hash(y))
        end
    end
end
