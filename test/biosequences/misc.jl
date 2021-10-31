@testset "Convert to String or Vector" begin
    seq = SimpleSeq("ACGUAAUUUCA")

    @test String(seq) == "ACGUAAUUUCA"

    seq = SimpleSeq(RNA[])
    @test isempty(seq)
end

@testset "Counting" begin
    for i in 1:3
        seq = random_simple(100)
        str = String(seq)
        @test count(isequal(RNA_A), seq) == count(isequal('A'), str)
        @test count(isambiguous, seq) == 0
        @test count(iscertain, seq) == length(seq)
        @test count(x -> x in (RNA_A, RNA_C), seq) == count(x -> x in "AC", str)

        @test isapprox(gc_content(seq), count(x -> x in "GC", str) / length(seq))
        @test n_ambiguous(seq) == 0
        @test n_certain(seq) == length(seq)
        @test n_gaps(seq) == 0
    end
end

@testset "Matches/mismatches" begin
    seq1 = random_simple(1000)
    seq2 = random_simple(1000)
    n_matches = sum(i == j for (i, j) in zip(seq1, seq2))
    @test mismatches(seq1, seq2) == length(seq1) - n_matches
    @test matches(seq1, seq2) == n_matches
end

@testset "Finding" begin
    @testset "Findfirst" begin
        seq = SimpleSeq("AUGCUGAUGAC")
        seq2 = SimpleSeq("AUGAUGAUGAUGUACA")
        @test findfirst(isequal(RNA_U), seq) == 2
        @test findfirst(isequal(missing), seq) === nothing
        @test findfirst(x -> true, empty(seq)) === nothing
        @test findfirst(x -> x == RNA_C, seq2) == 15
    end

    @testset "Findlast" begin
        seq = SimpleSeq("AUGCUGAUGAC")
        seq2 = SimpleSeq("AUGAUGAUGAUGUACA")
        @test findlast(isequal(RNA_U), seq) == 8
        @test findlast(isequal(missing), seq) === nothing
        @test findlast(x -> true, empty(seq)) === nothing
        @test findlast(x -> x == RNA_C, seq2) == 15
    end

    @testset "Findnext / prev" begin
        seq = SimpleSeq("AUGAUGAUGAUGUACA") # 16 nt

        @test findnext(x -> true, seq, 17) === nothing
        @test findprev(x -> true, seq, 0) === nothing

        @test_throws BoundsError findnext(x -> true, seq, 0)
        @test_throws BoundsError findprev(x -> true, seq, 17)

        @test findnext(isequal(RNA_U), seq, 6) == 8
        @test findprev(isequal(RNA_U), seq, 6) == 5
    end
end

@testset "Equality" begin
    @test SimpleSeq("AUGC") == SimpleSeq("AUGC")
    @test SimpleSeq("AUGC") != SimpleSeq("AUG")
    @test SimpleSeq("AUGC") != SimpleSeq("AUGCA")
    @test SimpleSeq("A") != SimpleSeq("U")

    @test !isless(SimpleSeq("GAC"), SimpleSeq("CAC"))
    @test isless(SimpleSeq("UG"), SimpleSeq("UGA"))
    @test isless(SimpleSeq("AGCUUA"), SimpleSeq("AGCUUU"))

    # This is particular for the RNA alphabet of SimpleSeq, not generic
    for i in 1:5
        seq1, seq2 = random_simple(20), random_simple(20)
        @test isless(seq1, seq2) == isless(String(seq1), String(seq2))
    end
end

@testset "Repetitive" begin
    @test isrepetitive(SimpleSeq("CCCCCCCCC"))
    @test !isrepetitive(SimpleSeq("CCCCCCCCA"))
    @test isrepetitive(SimpleSeq(RNA[]))

    @test isrepetitive(SimpleSeq("GAUGUCGAAAC"), 3)
    @test !isrepetitive(SimpleSeq("GAUGUCGAAAC"), 4)

    @test  isrepetitive(SimpleSeq(""))
    @test !isrepetitive(SimpleSeq(""), 1)
    @test  isrepetitive(SimpleSeq("A"))
    @test  isrepetitive(SimpleSeq("A"), 1)
    @test  isrepetitive(SimpleSeq("AAA"))
    @test !isrepetitive(SimpleSeq("ACGU"), 2)
    @test  isrepetitive(SimpleSeq("AAGU"), 2)
    @test  isrepetitive(SimpleSeq("ACCG"), 2)
    @test  isrepetitive(SimpleSeq("ACGG"), 2)
    @test !isrepetitive(SimpleSeq("ACGUCCGU"), 3)
    @test  isrepetitive(SimpleSeq("ACGCCCGU"), 3)

    @test  isrepetitive(SimpleSeq(""))
    @test !isrepetitive(SimpleSeq(""), 1)
    @test  isrepetitive(SimpleSeq("A"))
    @test  isrepetitive(SimpleSeq("A"), 1)
    @test  isrepetitive(SimpleSeq("AAA"))
    @test !isrepetitive(SimpleSeq("ACGU"), 2)
    @test  isrepetitive(SimpleSeq("AAGU"), 2)
    @test  isrepetitive(SimpleSeq("ACCG"), 2)
    @test  isrepetitive(SimpleSeq("ACGG"), 2)
    @test !isrepetitive(SimpleSeq("ACGUCCGU"), 3)
    @test  isrepetitive(SimpleSeq("ACGCCCGU"), 3)
end

@testset "Canonical" begin
    @test iscanonical(SimpleSeq("ACCG"))
    @test iscanonical(SimpleSeq("GCAC"))
    @test iscanonical(SimpleSeq("AAUU"))
    @test !iscanonical(SimpleSeq("UGGA"))
    @test !iscanonical(SimpleSeq("CGAU"))

    @test canonical(SimpleSeq("UGGA")) == SimpleSeq("UCCA")
    @test canonical(SimpleSeq("GCAC")) == SimpleSeq("GCAC")
    seq = SimpleSeq("CGAU")
    canonical!(seq)
    @test seq == SimpleSeq("AUCG")
end

@testset "Ispalindromic" begin
    @test  ispalindromic(SimpleSeq(""))
    @test !ispalindromic(SimpleSeq("A"))
    @test !ispalindromic(SimpleSeq("C"))
    @test  ispalindromic(SimpleSeq("AU"))
    @test  ispalindromic(SimpleSeq("CG"))
    @test !ispalindromic(SimpleSeq("AC"))
    @test !ispalindromic(SimpleSeq("UU"))
    @test  ispalindromic(SimpleSeq("ACGU"))

    @test  ispalindromic(SimpleSeq(""))
    @test !ispalindromic(SimpleSeq("A"))
    @test !ispalindromic(SimpleSeq("C"))
    @test  ispalindromic(SimpleSeq("AU"))
    @test  ispalindromic(SimpleSeq("CG"))
    @test !ispalindromic(SimpleSeq("AC"))
    @test !ispalindromic(SimpleSeq("UU"))
    @test  ispalindromic(SimpleSeq("ACGU"))
end

@testset "Has ambiguity" begin
    @test !hasambiguity(SimpleSeq("UAGUCGUGAG"))

    @test !hasambiguity(SimpleSeq(""))
    @test !hasambiguity(SimpleSeq("A"))
    @test !hasambiguity(SimpleSeq("ACGU"))

    @test !hasambiguity(SimpleSeq(""))
    @test !hasambiguity(SimpleSeq("A"))
    @test !hasambiguity(SimpleSeq("ACGU"))
end

@testset "Shuffle" begin
    function test_same(a, b)
        @test all(symbols(Alphabet(a))) do i
            count(isequal(i), a) == count(isequal(i), b)
        end
    end
    seq = SimpleSeq([RNA(i) for i in "AGCGUUAUGCUGAUUAGGAC"])
    seq2 = Random.shuffle(seq)
    test_same(seq, seq2)
    Random.shuffle!(seq)
    test_same(seq, seq2)
end

@testset "Reverse-complement" begin
    seq = SimpleSeq([RNA(i) for i in "UAGUUC"])
    @test reverse(seq) == SimpleSeq([RNA(i) for i in "CUUGAU"])
    @test complement(seq) == SimpleSeq([RNA(i) for i in "AUCAAG"])
    @test reverse_complement(seq) == reverse(complement(seq))

    reverse!(seq)
    @test seq == SimpleSeq([RNA(i) for i in "CUUGAU"])
    complement!(seq)
    @test seq == SimpleSeq([RNA(i) for i in "GAACUA"])
end

@testset "Ungap" begin
    seq = SimpleSeq([RNA(i) for i in "UAGUUC"])
    @test ungap(seq) == seq
    cp = copy(seq)
    @test ungap!(seq) == cp
end