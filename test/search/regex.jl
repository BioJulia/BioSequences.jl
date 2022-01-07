@testset "Regular Expression" begin
    # First test the biore syntax works as intended
    @test biore"A+[CG]N"d == BioRegex{DNA}("A+[CG]N")
    @test biore"A+[CG]N"d != BioRegex{DNA}("A*[CG]N")
    @test biore"MV?(I|L)W{3,5}$"a == BioRegex{AminoAcid}("MV?(I|L)W{3,5}\$")
    @test biore"MV?(I|L)W{3,5}$"a != BioRegex{AminoAcid}("MV?(I|L)W{3,4}\$")
    @test biore"A+UUC*"r == BioRegex{RNA}("A+UUC*")
    @test biore"A+UC*"r != BioRegex{RNA}("A+UUC*")

    @test isa(biore"A+"d, BioSequences.RE.Regex{DNA})
    @test isa(biore"A+"r, BioSequences.RE.Regex{RNA})
    @test isa(biore"A+"a, BioSequences.RE.Regex{AminoAcid})
    @test isa(biore"A+"dna, BioSequences.RE.Regex{DNA})
    @test isa(biore"A+"rna, BioSequences.RE.Regex{RNA})
    @test isa(biore"A+"aa, BioSequences.RE.Regex{AminoAcid})
    @test_throws Exception eval(:(biore"A+"))
    @test_throws Exception eval(:(biore"A+"foo))
    @test string(biore"A+"dna) == "biore\"A+\"dna"
    @test string(biore"A+"rna) == "biore\"A+\"rna"
    @test string(biore"A+"aa) == "biore\"A+\"aa"

    @test  occursin(biore"A"d, dna"A")
    @test !occursin(biore"A"d, dna"C")
    @test !occursin(biore"A"d, dna"G")
    @test !occursin(biore"A"d, dna"T")
    @test  occursin(biore"N"d, dna"A")
    @test  occursin(biore"N"d, dna"C")
    @test  occursin(biore"N"d, dna"G")
    @test  occursin(biore"N"d, dna"T")
    @test  occursin(biore"[AT]"d, dna"A")
    @test !occursin(biore"[AT]"d, dna"C")
    @test !occursin(biore"[AT]"d, dna"G")
    @test  occursin(biore"[AT]"d, dna"T")
    @test !occursin(biore"[^AT]"d, dna"A")
    @test  occursin(biore"[^AT]"d, dna"C")
    @test  occursin(biore"[^AT]"d, dna"G")
    @test !occursin(biore"[^AT]"d, dna"T")

    re = biore"^A(C+G*)(T{2,})N$"d
    @test !occursin(re, dna"AC")
    @test !occursin(re, dna"AGTT")
    @test !occursin(re, dna"CCGTT")
    @test !occursin(re, dna"ACTT")
    @test !occursin(re, dna"ACTTGT")
    @test  occursin(re, dna"ACGTTA")
    @test  occursin(re, dna"ACGTTT")
    @test  occursin(re, dna"ACCGGTTT")
    @test  occursin(re, dna"ACCGGTTT")
    @test  occursin(re, dna"ACCGGTTTA")
    @test  occursin(re, dna"ACCGGTTTG")

    @test matched(match(re, dna"ACCGTTTTA")) == dna"ACCGTTTTA"
    @test captured(match(re, dna"ACCGTTTTA"))[1] == dna"CCG"
    @test captured(match(re, dna"ACCGTTTTA"))[2] == dna"TTTT"

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
    @test findfirst(biore"A+"d, dna"ACGTAAT") == 1:1
    @test findfirst(biore"A+"d, dna"ACGTAAT", 1) == 1:1
    @test findfirst(biore"A+"d, dna"ACGTAAT", 2) == 5:6
    @test findfirst(biore"A+"d, dna"ACGTAAT", 7) === nothing

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
    matchall(pat, seq, overlap=true) =
        collect(map(matched, eachmatch(pat, seq, overlap)))

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
    @test  occursin(biore"U(A[AG]|GA)$"r, rna"AUUGUAUGA")
    @test !occursin(biore"U(A[AG]|GA)$"r, rna"AUUGUAUGG")
    @test  occursin(biore"T+[NQ]A?P"a, aa"MTTQAPMFTQPL")
    @test  occursin(biore"T+[NQ]A?P"a, aa"MTTAAPMFTQPL")
    @test !occursin(biore"T+[NQ]A?P"a, aa"MTTAAPMFSQPL")

    # PROSITE
    @test  occursin(prosite"[AC]-x-V-x(4)-{ED}", aa"ADVAARRK")
    @test  occursin(prosite"[AC]-x-V-x(4)-{ED}", aa"CPVAARRK")
    @test !occursin(prosite"[AC]-x-V-x(4)-{ED}", aa"ADVAARRE")
    @test !occursin(prosite"[AC]-x-V-x(4)-{ED}", aa"CPVAARK")
    @test  occursin(prosite"<[AC]-x-V-x(4)-{ED}>", aa"ADVAARRK")
    @test !occursin(prosite"<[AC]-x-V-x(4)-{ED}>", aa"AADVAARRK")
    @test !occursin(prosite"<[AC]-x-V-x(4)-{ED}>", aa"ADVAARRKA")
end
