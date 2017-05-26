@testset "Regular Expression" begin
    @test isa(biore"A+"d, Seq.RE.Regex{DNA})
    @test isa(biore"A+"r, Seq.RE.Regex{RNA})
    @test isa(biore"A+"a, Seq.RE.Regex{AminoAcid})
    @test isa(biore"A+"dna, Seq.RE.Regex{DNA})
    @test isa(biore"A+"rna, Seq.RE.Regex{RNA})
    @test isa(biore"A+"aa, Seq.RE.Regex{AminoAcid})
    @test_throws Exception eval(:(biore"A+"))
    @test_throws Exception eval(:(biore"A+"foo))
    @test string(biore"A+"dna) == "biore\"A+\"dna"
    @test string(biore"A+"rna) == "biore\"A+\"rna"
    @test string(biore"A+"aa) == "biore\"A+\"aa"

    @test  ismatch(biore"A"d, dna"A")
    @test !ismatch(biore"A"d, dna"C")
    @test !ismatch(biore"A"d, dna"G")
    @test !ismatch(biore"A"d, dna"T")
    @test  ismatch(biore"N"d, dna"A")
    @test  ismatch(biore"N"d, dna"C")
    @test  ismatch(biore"N"d, dna"G")
    @test  ismatch(biore"N"d, dna"T")
    @test  ismatch(biore"[AT]"d, dna"A")
    @test !ismatch(biore"[AT]"d, dna"C")
    @test !ismatch(biore"[AT]"d, dna"G")
    @test  ismatch(biore"[AT]"d, dna"T")
    @test !ismatch(biore"[^AT]"d, dna"A")
    @test  ismatch(biore"[^AT]"d, dna"C")
    @test  ismatch(biore"[^AT]"d, dna"G")
    @test !ismatch(biore"[^AT]"d, dna"T")

    re = biore"^A(C+G*)(T{2,})N$"d
    @test !ismatch(re, dna"AC")
    @test !ismatch(re, dna"AGTT")
    @test !ismatch(re, dna"CCGTT")
    @test !ismatch(re, dna"ACTT")
    @test !ismatch(re, dna"ACTTGT")
    @test  ismatch(re, dna"ACGTTA")
    @test  ismatch(re, dna"ACGTTT")
    @test  ismatch(re, dna"ACCGGTTT")
    @test  ismatch(re, dna"ACCGGTTT")
    @test  ismatch(re, dna"ACCGGTTTA")
    @test  ismatch(re, dna"ACCGGTTTG")

    @test matched(match(re, dna"ACCGTTTTA")) == dna"ACCGTTTTA"
    @test get(captured(match(re, dna"ACCGTTTTA"))[1]) == dna"CCG"
    @test get(captured(match(re, dna"ACCGTTTTA"))[2]) == dna"TTTT"

    # greedy
    @test matched(match(biore"A*"d, dna"AAA")) == dna"AAA"
    @test matched(match(biore"A+"d, dna"AAA")) == dna"AAA"
    @test matched(match(biore"A?"d, dna"AAA")) == dna"A"
    @test matched(match(biore"A{2}"d, dna"AAA")) == dna"AA"
    @test matched(match(biore"A{2,}"d, dna"AAA")) == dna"AAA"
    @test matched(match(biore"A{2,4}"d, dna"AAA")) == dna"AAA"

    # lazy
    @test matched(match(biore"A*?"d, dna"AAA")) == dna""
    @test matched(match(biore"A+?"d, dna"AAA")) == dna"A"
    @test matched(match(biore"A??"d, dna"AAA")) == dna""
    @test matched(match(biore"A{2}?"d, dna"AAA")) == dna"AA"
    @test matched(match(biore"A{2,}?"d, dna"AAA")) == dna"AA"
    @test matched(match(biore"A{2,4}?"d, dna"AAA")) == dna"AA"

    # search
    @test search(dna"ACGTAAT", biore"A+"d) == 1:1
    @test search(dna"ACGTAAT", biore"A+"d, 1) == 1:1
    @test search(dna"ACGTAAT", biore"A+"d, 2) == 5:6
    @test search(dna"ACGTAAT", biore"A+"d, 7) == 0:-1

    # eachmatch
    matches = [dna"CG", dna"GC", dna"GC", dna"CG"]
    for (i, m) in enumerate(eachmatch(biore"GC|CG"d, dna"ACGTTATGCATGGCG"))
        @test matched(m) == matches[i]
    end
    matches = [dna"CG", dna"GC", dna"GC"]
    for (i, m) in enumerate(collect(eachmatch(biore"GC|CG"d, dna"ACGTTATGCATGGCG", false)))
        @test matched(m) == matches[i]
    end

    # matchall
    @test matchall(biore"A*"d, dna"") == [dna""]
    @test matchall(biore"A*"d, dna"AAA") == [
        dna"AAA", dna"AA", dna"A", dna"",
        dna"AA",  dna"A",  dna"",
        dna"A",   dna""]
    @test matchall(biore"AC*G*T"d, dna"ACCGGGT") == [dna"ACCGGGT"]

    @test matchall(biore"A*"d, dna"", false) == [dna""]
    @test matchall(biore"A*"d, dna"AAA", false) == [dna"AAA"]
    @test matchall(biore"AC*G*T"d, dna"ACCGGGT", false) == [dna"ACCGGGT"]

    # RNA and Amino acid
    @test  ismatch(biore"U(A[AG]|GA)$"r, rna"AUUGUAUGA")
    @test !ismatch(biore"U(A[AG]|GA)$"r, rna"AUUGUAUGG")
    @test  ismatch(biore"T+[NQ]A?P"a, aa"MTTQAPMFTQPL")
    @test  ismatch(biore"T+[NQ]A?P"a, aa"MTTAAPMFTQPL")
    @test !ismatch(biore"T+[NQ]A?P"a, aa"MTTAAPMFSQPL")

    # PROSITE
    @test  ismatch(prosite"[AC]-x-V-x(4)-{ED}", aa"ADVAARRK")
    @test  ismatch(prosite"[AC]-x-V-x(4)-{ED}", aa"CPVAARRK")
    @test !ismatch(prosite"[AC]-x-V-x(4)-{ED}", aa"ADVAARRE")
    @test !ismatch(prosite"[AC]-x-V-x(4)-{ED}", aa"CPVAARK")
    @test  ismatch(prosite"<[AC]-x-V-x(4)-{ED}>", aa"ADVAARRK")
    @test !ismatch(prosite"<[AC]-x-V-x(4)-{ED}>", aa"AADVAARRK")
    @test !ismatch(prosite"<[AC]-x-V-x(4)-{ED}>", aa"ADVAARRKA")
end
