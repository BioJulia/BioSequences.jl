@testset "Comparisons" begin
    @testset "Equality" begin
        function check_seq_kmer_equality(len)
            a = DNAMer(random_dna_kmer(len))
            b = LongDNASeq(a)
            return a == b && b == a
        end
        
        function check_seq_bigkmer_equality(len)
            a = BigDNAMer(random_dna_kmer(len))
            b = LongDNASeq(a)
            return a == b && b == a
        end

        for len in [1, 10, 32, 64]
            if len <= 32
                @test all(Bool[check_seq_kmer_equality(len) for _ in 1:reps])
            end
            @test all(Bool[check_seq_bigkmer_equality(len) for _ in 1:reps])
        end

        # True negatives
        @test DNAMer("ACG") != RNAMer("ACG")
        @test DNAMer("T")   != RNAMer("U")
        @test DNAMer("AC")  != DNAMer("AG")
        @test RNAMer("AC")  != RNAMer("AG")
        
        @test BigDNAMer("ACG") != BigRNAMer("ACG")
        @test BigDNAMer("T")   != BigRNAMer("U")
        @test BigDNAMer("AC")  != BigDNAMer("AG")
        @test BigRNAMer("AC")  != BigRNAMer("AG")

        @test DNAMer("ACG") != rna"ACG"
        @test DNAMer("T")   != rna"U"
        @test DNAMer("AC")  != dna"AG"
        @test RNAMer("AC")  != rna"AG"
        
        @test BigDNAMer("ACG") != rna"ACG"
        @test BigDNAMer("T")   != rna"U"
        @test BigDNAMer("AC")  != dna"AG"
        @test BigRNAMer("AC")  != rna"AG"

        @test rna"ACG" != DNAMer("ACG")
        @test rna"U"   != DNAMer("T")
        @test dna"AG"  != DNAMer("AC")
        @test rna"AG"  != RNAMer("AC")
        
        @test rna"ACG" != BigDNAMer("ACG")
        @test rna"U"   != BigDNAMer("T")
        @test dna"AG"  != BigDNAMer("AC")
        @test rna"AG"  != BigRNAMer("AC")
    end

    @testset "Inequality" begin
        for len in [1, 10, 32, 64]
            if len <= 32
                @test  isless(DNAMer{1}(UInt64(0)), DNAMer{1}(UInt64(1)))
                @test !isless(DNAMer{1}(UInt64(0)), DNAMer{1}(UInt64(0)))
                @test !isless(DNAMer{1}(UInt64(1)), DNAMer{1}(UInt64(0)))
                
                @test  isless(RNAMer{1}(UInt64(0)), RNAMer{1}(UInt64(1)))
                @test !isless(RNAMer{1}(UInt64(0)), RNAMer{1}(UInt64(0)))
                @test !isless(RNAMer{1}(UInt64(1)), RNAMer{1}(UInt64(0)))
            end
            
            @test  isless(BigDNAMer{1}(UInt128(0)), BigDNAMer{1}(UInt128(1)))
            @test !isless(BigDNAMer{1}(UInt128(0)), BigDNAMer{1}(UInt128(0)))
            @test !isless(BigDNAMer{1}(UInt128(1)), BigDNAMer{1}(UInt128(0)))
            
            @test  isless(BigRNAMer{1}(UInt128(0)), BigRNAMer{1}(UInt128(1)))
            @test !isless(BigRNAMer{1}(UInt128(0)), BigRNAMer{1}(UInt128(0)))
            @test !isless(BigRNAMer{1}(UInt128(1)), BigRNAMer{1}(UInt128(0)))
        end
    end

    @testset "Hash" begin
        kmers = map(DNAMer, ["AAAA", "AACT", "ACGT", "TGCA"])
        for x in kmers, y in kmers
            @test (x == y) == (hash(x) == hash(y))
        end
        
        bigkmers = map(BigDNAMer, ["AAAA", "AACT", "ACGT", "TGCA"])
        for x in bigkmers, y in bigkmers
            @test (x == y) == (hash(x) == hash(y))
        end
        
        kmers = map(RNAMer, ["AAAA", "AACU", "ACGU", "UGCA"])
        for x in kmers, y in kmers
            @test (x == y) == (hash(x) == hash(y))
        end
        
        bigkmers = map(BigRNAMer, ["AAAA", "AACU", "ACGU", "UGCA"])
        for x in bigkmers, y in bigkmers
            @test (x == y) == (hash(x) == hash(y))
        end
    end
end
