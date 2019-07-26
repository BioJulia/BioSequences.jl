@testset "Access and Iterations" begin
    dna_kmer = mer"ACTG"dna
    rna_kmer = mer"ACUG"rna
    big_dna_kmer = bigmer"ACTG"dna
    big_rna_kmer = bigmer"ACUG"rna

    @testset "Access DNA Kmer" begin
        @test dna_kmer[1] == DNA_A
        @test dna_kmer[2] == DNA_C
        @test dna_kmer[3] == DNA_T
        @test dna_kmer[4] == DNA_G
        
        @test big_dna_kmer[1] == DNA_A
        @test big_dna_kmer[2] == DNA_C
        @test big_dna_kmer[3] == DNA_T
        @test big_dna_kmer[4] == DNA_G

        # Access indexes out of bounds
        @test_throws Exception dna_kmer[-1]
        @test_throws Exception dna_kmer[0]
        @test_throws Exception dna_kmer[5]
        @test_throws Exception getindex(dna_kmer,-1)
        @test_throws Exception getindex(dna_kmer, 0)
        @test_throws Exception getindex(dna_kmer, 5)
        
        @test_throws Exception big_dna_kmer[-1]
        @test_throws Exception big_dna_kmer[0]
        @test_throws Exception big_dna_kmer[5]
        @test_throws Exception getindex(big_dna_kmer,-1)
        @test_throws Exception getindex(big_dna_kmer, 0)
        @test_throws Exception getindex(big_dna_kmer, 5)
    end

    @testset "Iteration through DNA Kmer" begin
        @test iterate(DNAMer("ACTG")) == (DNA_A, 2)
        @test iterate(BigDNAMer("ACTG")) == (DNA_A, 2)

        @test iterate(DNAMer("ACTG"), 1) == (DNA_A, 2)
        @test iterate(DNAMer("ACTG"), 4) == (DNA_G, 5)
        @test iterate(BigDNAMer("ACTG"), 1) == (DNA_A, 2)
        @test iterate(BigDNAMer("ACTG"), 4) == (DNA_G, 5)

        @test iterate(DNAMer("ACTG"), 1)  !== nothing
        @test iterate(DNAMer("ACTG"), 4)  !== nothing
        @test iterate(DNAMer("ACTG"), 5)  === nothing
        @test iterate(DNAMer("ACTG"), -1) !== nothing
        
        @test iterate(BigDNAMer("ACTG"), 1)  !== nothing
        @test iterate(BigDNAMer("ACTG"), 4)  !== nothing
        @test iterate(BigDNAMer("ACTG"), 5)  === nothing
        @test iterate(BigDNAMer("ACTG"), -1) !== nothing

        dna_vec = [DNA_A, DNA_C, DNA_T, DNA_G]
        @test all([nt === dna_vec[i] for (i, nt) in enumerate(dna_kmer)])
        @test all([nt === dna_vec[i] for (i, nt) in enumerate(big_dna_kmer)])
    end

    @testset "Access RNA Kmer" begin
        @test rna_kmer[1] == RNA_A
        @test rna_kmer[2] == RNA_C
        @test rna_kmer[3] == RNA_U
        @test rna_kmer[4] == RNA_G
        
        @test big_rna_kmer[1] == RNA_A
        @test big_rna_kmer[2] == RNA_C
        @test big_rna_kmer[3] == RNA_U
        @test big_rna_kmer[4] == RNA_G

        # Access indexes out of bounds
        @test_throws Exception rna_kmer[-1]
        @test_throws Exception rna_kmer[0]
        @test_throws Exception rna_kmer[5]
        @test_throws Exception getindex(rna_kmer, -1)
        @test_throws Exception getindex(rna_kmer, 0)
        @test_throws Exception getindex(rna_kmer, 5)
        
        @test_throws Exception big_rna_kmer[-1]
        @test_throws Exception big_rna_kmer[0]
        @test_throws Exception big_rna_kmer[5]
        @test_throws Exception getindex(big_rna_kmer, -1)
        @test_throws Exception getindex(big_rna_kmer, 0)
        @test_throws Exception getindex(big_rna_kmer, 5)
    end

    @testset "Iteration through RNA Kmer" begin
        @test iterate(RNAMer("ACUG")) == (RNA_A, 2)
        @test iterate(BigRNAMer("ACUG")) == (RNA_A, 2)

        @test iterate(RNAMer("ACUG"), 1) == (RNA_A, 2)
        @test iterate(RNAMer("ACUG"), 4) == (RNA_G, 5)
        @test iterate(BigRNAMer("ACUG"), 1) == (RNA_A, 2)
        @test iterate(BigRNAMer("ACUG"), 4) == (RNA_G, 5)

        @test iterate(RNAMer("ACUG"), 1)  !== nothing
        @test iterate(RNAMer("ACUG"), 4)  !== nothing
        @test iterate(RNAMer("ACUG"), 5)  === nothing
        @test iterate(RNAMer("ACUG"), -1) !== nothing
        
        @test iterate(BigRNAMer("ACUG"), 1)  !== nothing
        @test iterate(BigRNAMer("ACUG"), 4)  !== nothing
        @test iterate(BigRNAMer("ACUG"), 5)  === nothing
        @test iterate(BigRNAMer("ACUG"), -1) !== nothing

        rna_vec = [RNA_A, RNA_C, RNA_U, RNA_G]
        @test all([nt === rna_vec[i] for (i, nt) in enumerate(rna_kmer)])
        @test all([nt === rna_vec[i] for (i, nt) in enumerate(big_rna_kmer)])
    end
end
