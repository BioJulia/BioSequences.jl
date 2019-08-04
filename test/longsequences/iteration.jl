@testset "Iteration" begin
    dna_seq = dna"ACTG"
    dna_vec = [DNA_A, DNA_C, DNA_T, DNA_G]
    @test all([nt === dna_vec[i] for (i, nt) in enumerate(dna_seq)])

    rna_seq = rna"ACUG"
    rna_vec = [RNA_A, RNA_C, RNA_U, RNA_G]
    @test all([nt === rna_vec[i] for (i, nt) in enumerate(rna_seq)])

    aa_seq = aa"ARNPS"
    aa_vec = [AA_A, AA_R, AA_N, AA_P, AA_S]
    @test all([aa == aa_vec[i] for (i, aa) in enumerate(aa_seq)])
end
