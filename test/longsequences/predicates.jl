@testset "Predicates" begin
    # ispalindromic
    @test  ispalindromic(dna"")
    @test !ispalindromic(dna"A")
    @test !ispalindromic(dna"C")
    @test  ispalindromic(dna"AT")
    @test  ispalindromic(dna"CG")
    @test !ispalindromic(dna"AC")
    @test !ispalindromic(dna"TT")
    @test  ispalindromic(dna"ANT")
    @test  ispalindromic(dna"ACGT")
    @test !ispalindromic(dna"ACNT")

    @test  ispalindromic(rna"")
    @test !ispalindromic(rna"A")
    @test !ispalindromic(rna"C")
    @test  ispalindromic(rna"AU")
    @test  ispalindromic(rna"CG")
    @test !ispalindromic(rna"AC")
    @test !ispalindromic(rna"UU")
    @test  ispalindromic(rna"ANU")
    @test  ispalindromic(rna"ACGU")
    @test !ispalindromic(rna"ACNU")

    @test_throws Exception ispalindromic(aa"PQ")

    # hasambiguity
    @test !hasambiguity(dna"")
    @test !hasambiguity(dna"A")
    @test  hasambiguity(dna"N")
    @test !hasambiguity(dna"ACGT")
    @test  hasambiguity(dna"ANGT")

    @test !hasambiguity(rna"")
    @test !hasambiguity(rna"A")
    @test  hasambiguity(rna"N")
    @test !hasambiguity(rna"ACGU")
    @test  hasambiguity(rna"ANGU")

    @test !hasambiguity(aa"")
    @test !hasambiguity(aa"A")
    @test !hasambiguity(aa"P")
    @test  hasambiguity(aa"B")
    @test  hasambiguity(aa"X")
    @test !hasambiguity(aa"ARNDCQEGHILKMFPSTWYVOU")
    @test  hasambiguity(aa"ARXDCQEGHILKMFPSTWYVOU")

    # isrepetitive
    @test  isrepetitive(dna"")
    @test !isrepetitive(dna"", 1)
    @test  isrepetitive(dna"A")
    @test  isrepetitive(dna"A", 1)
    @test  isrepetitive(dna"AAA")
    @test !isrepetitive(dna"ACGT", 2)
    @test  isrepetitive(dna"AAGT", 2)
    @test  isrepetitive(dna"ACCG", 2)
    @test  isrepetitive(dna"ACGG", 2)
    @test !isrepetitive(dna"ACGTCCGT", 3)
    @test  isrepetitive(dna"ACGCCCGT", 3)

    @test  isrepetitive(rna"")
    @test !isrepetitive(rna"", 1)
    @test  isrepetitive(rna"A")
    @test  isrepetitive(rna"A", 1)
    @test  isrepetitive(rna"AAA")
    @test !isrepetitive(rna"ACGU", 2)
    @test  isrepetitive(rna"AAGU", 2)
    @test  isrepetitive(rna"ACCG", 2)
    @test  isrepetitive(rna"ACGG", 2)
    @test !isrepetitive(rna"ACGUCCGU", 3)
    @test  isrepetitive(rna"ACGCCCGU", 3)

    @test  isrepetitive(aa"")
    @test !isrepetitive(aa"PGQQ")
    @test  isrepetitive(aa"PGQQ", 2)
    @test !isrepetitive(aa"PPQQ", 3)
    @test  isrepetitive(aa"PPPQQ", 3)

    # iscanonical
    @test iscanonical(dna"TCA")
    @test !iscanonical(dna"TGA")

    @test !hasprematurestop(dna"ATGCTTAAACCTTTGACGTAG")
    @test hasprematurestop(dna"ATGCTTAAACCTTTGACGTAGTAG") # stop (TAG) before last codon
    @test hasprematurestop(dna"TAGATGCTTAAACCTTTGACGTAG") # stop (TAG) at the beginning
    
end
