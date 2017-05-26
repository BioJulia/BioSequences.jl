@testset "Access and Iterations" begin
    dna_kmer = DNAKmer("ACTG")
    rna_kmer = RNAKmer("ACUG")

    @testset "Access DNA Kmer" begin
        @test dna_kmer[1] == DNA_A
        @test dna_kmer[2] == DNA_C
        @test dna_kmer[3] == DNA_T
        @test dna_kmer[4] == DNA_G

        # Access indexes out of bounds
        @test_throws Exception dna_kmer[-1]
        @test_throws Exception dna_kmer[0]
        @test_throws Exception dna_kmer[5]
        @test_throws Exception getindex(dna_kmer,-1)
        @test_throws Exception getindex(dna_kmer, 0)
        @test_throws Exception getindex(dna_kmer, 5)
    end

    @testset "Iteration through DNA Kmer" begin
        @test start(DNAKmer("ACTG")) == 1
        @test start(DNAKmer(""))     == 1

        @test next(DNAKmer("ACTG"), 1) == (DNA_A, 2)
        @test next(DNAKmer("ACTG"), 4) == (DNA_G, 5)

        @test  done(DNAKmer(""), 1)
        @test !done(DNAKmer("ACTG"), 1)
        @test !done(DNAKmer("ACTG"), 4)
        @test  done(DNAKmer("ACTG"), 5)
        @test !done(DNAKmer("ACTG"), -1)


        dna_vec = [DNA_A, DNA_C, DNA_T, DNA_G]
        @test all([nt === dna_vec[i] for (i, nt) in enumerate(dna_kmer)])
    end

    @testset "Access RNA Kmer" begin
        @test rna_kmer[1] == RNA_A
        @test rna_kmer[2] == RNA_C
        @test rna_kmer[3] == RNA_U
        @test rna_kmer[4] == RNA_G

        # Access indexes out of bounds
        @test_throws Exception rna_kmer[-1]
        @test_throws Exception rna_kmer[0]
        @test_throws Exception rna_kmer[5]
        @test_throws Exception getindex(rna_kmer, -1)
        @test_throws Exception getindex(rna_kmer, 0)
        @test_throws Exception getindex(rna_kmer, 5)
    end

    @testset "Iteration through RNA Kmer" begin
        @test start(RNAKmer("ACUG")) == 1
        @test start(RNAKmer(""))     == 1

        @test next(RNAKmer("ACUG"), 1) == (RNA_A, 2)
        @test next(RNAKmer("ACUG"), 4) == (RNA_G, 5)

        @test  done(RNAKmer(""), 1)
        @test !done(RNAKmer("ACUG"), 1)
        @test !done(RNAKmer("ACUG"), 4)
        @test  done(RNAKmer("ACUG"), 5)
        @test !done(RNAKmer("ACUG"), -1)

        rna_vec = [RNA_A, RNA_C, RNA_U, RNA_G]
        @test all([nt === rna_vec[i] for (i, nt) in enumerate(rna_kmer)])
    end
end
