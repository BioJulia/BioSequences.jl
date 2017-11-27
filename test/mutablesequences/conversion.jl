@testset "Constructing empty sequences" begin
    @test DNASequence() == MutableBioSequence(DNA)
    @test RNASequence() == MutableBioSequence(RNA)
    @test AminoAcidSequence() == MutableBioSequence(AminoAcid)
    @test CharSequence() == MutableBioSequence(Char)
end

@testset "Constructing uninitialized sequences" begin
    @test isa(MutableBioSequence{DNAAlphabet{2}}(0), MutableBioSequence)
    @test isa(MutableBioSequence{DNAAlphabet{4}}(10), MutableBioSequence)
    @test isa(MutableBioSequence{RNAAlphabet{2}}(0), MutableBioSequence)
    @test isa(MutableBioSequence{RNAAlphabet{4}}(10), MutableBioSequence)
    @test isa(MutableBioSequence{AminoAcidAlphabet}(10), MutableBioSequence)
end

@testset "Conversion from/to strings" begin
    # Check that sequences in strings survive round trip conversion:
    #   String → MutableBioSequence → String
    function test_string_construction(A::Type, seq::AbstractString)
        @test convert(String, MutableBioSequence{A}(seq)) == uppercase(seq)
    end

    function test_string_parse(A::Type, seq::AbstractString)
        @test parse(MutableBioSequence{A}, seq) == MutableBioSequence{A}(seq)
    end

    for len in [0, 1, 2, 3, 10, 32, 1000, 10000]
        test_string_construction(DNAAlphabet{4}, random_dna(len))
        test_string_construction(DNAAlphabet{4}, lowercase(random_dna(len)))
        test_string_construction(RNAAlphabet{4}, lowercase(random_rna(len)))
        test_string_construction(RNAAlphabet{4}, random_rna(len))
        test_string_construction(AminoAcidAlphabet, random_aa(len))
        test_string_construction(AminoAcidAlphabet, lowercase(random_aa(len)))

        test_string_parse(DNAAlphabet{4}, random_dna(len))
        test_string_parse(DNAAlphabet{4}, lowercase(random_dna(len)))
        test_string_parse(RNAAlphabet{4}, lowercase(random_rna(len)))
        test_string_parse(RNAAlphabet{4}, random_rna(len))
        test_string_parse(AminoAcidAlphabet, random_aa(len))
        test_string_parse(AminoAcidAlphabet, lowercase(random_aa(len)))

        probs = [0.25, 0.25, 0.25, 0.25, 0.00]
        test_string_construction(DNAAlphabet{2}, random_dna(len, probs))
        test_string_construction(DNAAlphabet{2}, lowercase(random_dna(len, probs)))
        test_string_construction(RNAAlphabet{2}, random_rna(len, probs))
        test_string_construction(RNAAlphabet{2}, lowercase(random_rna(len, probs)))

        test_string_parse(DNAAlphabet{2}, random_dna(len, probs))
        test_string_parse(DNAAlphabet{2}, lowercase(random_dna(len, probs)))
        test_string_parse(RNAAlphabet{2}, random_rna(len, probs))
        test_string_parse(RNAAlphabet{2}, lowercase(random_rna(len, probs)))
    end

    # non-standard string literal
    @test isa(dna"ACGTMRWSYKVHDBN-", DNASequence)
    @test isa(rna"ACGUMRWSYKVHDBN-", RNASequence)
    @test isa(aa"ARNDCQEGHILKMFPSTWYVBJZXOU*-", AminoAcidSequence)
    @test isa(char"いろは αβγ 甲乙丙", CharSequence)

    # Non-nucleotide characters should throw
    @test_throws Exception DNASequence("ACCNNCATTTTTTAGATXATAG")
    @test_throws Exception RNASequence("ACCNNCATTTTTTAGATXATAG")
    @test_throws Exception AminoAcidSequence("ATGHLMY@ZACAGNM")
end

@testset "Construction from vectors" begin
    function test_vector_construction(A, seq::AbstractString)
        T = eltype(A)
        xs = T[convert(T, c) for c in seq]
        @test MutableBioSequence{A}(xs) == MutableBioSequence{A}(seq)
    end

    for len in [0, 1, 10, 32, 1000, 10000]
        test_vector_construction(DNAAlphabet{4}, random_dna(len))
        test_vector_construction(RNAAlphabet{4}, random_rna(len))
        test_vector_construction(AminoAcidAlphabet, random_aa(len))

        probs = [0.25, 0.25, 0.25, 0.25, 0.00]
        test_vector_construction(DNAAlphabet{2}, random_dna(len, probs))
        test_vector_construction(RNAAlphabet{2}, random_rna(len, probs))
    end
end

@testset "Conversion between 2-bit and 4-bit encodings" begin
    function test_conversion(A1, A2, seq)
        @test convert(MutableBioSequence{A1}, MutableBioSequence{A2}(seq)) == convert(MutableBioSequence{A1}, seq)
    end

    test_conversion(DNAAlphabet{2}, DNAAlphabet{4}, "")
    test_conversion(DNAAlphabet{4}, DNAAlphabet{2}, "")
    test_conversion(RNAAlphabet{4}, RNAAlphabet{2}, "")
    test_conversion(RNAAlphabet{2}, RNAAlphabet{4}, "")

    test_conversion(DNAAlphabet{2}, DNAAlphabet{4}, "ACGT")
    test_conversion(DNAAlphabet{4}, DNAAlphabet{2}, "ACGT")
    test_conversion(RNAAlphabet{4}, RNAAlphabet{2}, "ACGU")
    test_conversion(RNAAlphabet{2}, RNAAlphabet{4}, "ACGU")

    test_conversion(DNAAlphabet{2}, DNAAlphabet{4}, "ACGT"^100)
    test_conversion(DNAAlphabet{4}, DNAAlphabet{2}, "ACGT"^100)
    test_conversion(RNAAlphabet{4}, RNAAlphabet{2}, "ACGU"^100)
    test_conversion(RNAAlphabet{2}, RNAAlphabet{4}, "ACGU"^100)

    # ambiguous nucleotides cannot be stored in 2-bit encoding
    EncodeError = BioSequences.EncodeError
    @test_throws EncodeError convert(MutableBioSequence{DNAAlphabet{2}}, dna"AN")
    @test_throws EncodeError convert(MutableBioSequence{RNAAlphabet{2}}, rna"AN")

    # test promotion
    a = MutableBioSequence{DNAAlphabet{2}}("ATCG")
    b = MutableBioSequence{DNAAlphabet{4}}("ATCG")
    c = MutableBioSequence{RNAAlphabet{2}}("AUCG")
    d = MutableBioSequence{RNAAlphabet{4}}("AUCG")

    @test typeof(promote(a, b)) == Tuple{MutableBioSequence{DNAAlphabet{4}},MutableBioSequence{DNAAlphabet{4}}}
    @test typeof(promote(c, d)) == Tuple{MutableBioSequence{RNAAlphabet{4}},MutableBioSequence{RNAAlphabet{4}}}
    @test_throws ErrorException typeof(promote(a, d))
    @test_throws ErrorException typeof(promote(a, b, d))
end

@testset "Conversion between RNA and DNA" begin
    @test convert(RNASequence, DNASequence("ACGTN")) == rna"ACGUN"
    @test convert(DNASequence, RNASequence("ACGUN")) == dna"ACGTN"
end

@testset "Conversion to Matrices" begin
    dna = [dna"AAA", dna"TTT", dna"CCC", dna"GGG"]
    dnathrow = [dna"AAA", dna"TTTAAA", dna"CCC", dna"GGG"]

    rna = [rna"AAA", rna"UUU", rna"CCC", rna"GGG"]
    rnathrow = [rna"AAA", rna"UUU", rna"CCCUUU", rna"GGG"]

    prot = [aa"AMG", aa"AMG", aa"AMG", aa"AMG"]

    sitemajdna = [
        DNA_A  DNA_A  DNA_A
        DNA_T  DNA_T  DNA_T
        DNA_C  DNA_C  DNA_C
        DNA_G  DNA_G  DNA_G
    ]
    seqmajdna = [
        DNA_A  DNA_T  DNA_C  DNA_G
        DNA_A  DNA_T  DNA_C  DNA_G
        DNA_A  DNA_T  DNA_C  DNA_G
    ]
    sitemajrna = [
        RNA_A  RNA_A  RNA_A
        RNA_U  RNA_U  RNA_U
        RNA_C  RNA_C  RNA_C
        RNA_G  RNA_G  RNA_G
    ]
    seqmajrna = [
        RNA_A  RNA_U  RNA_C  RNA_G
        RNA_A  RNA_U  RNA_C  RNA_G
        RNA_A  RNA_U  RNA_C  RNA_G
    ]
    sitemajnucint = [
        0x01  0x01  0x01
        0x08  0x08  0x08
        0x02  0x02  0x02
        0x04  0x04  0x04
    ]
    seqmajnucint = [
        0x01  0x08  0x02  0x04
        0x01  0x08  0x02  0x04
        0x01  0x08  0x02  0x04
    ]
    sitemajaa = [
        AA_A AA_M AA_G
        AA_A AA_M AA_G
        AA_A AA_M AA_G
        AA_A AA_M AA_G
    ]
    seqmajaa = [
        AA_A AA_A AA_A AA_A
        AA_M AA_M AA_M AA_M
        AA_G AA_G AA_G AA_G
    ]

    @test seqmatrix(dna, :site) == sitemajdna
    @test seqmatrix(rna, :site) == sitemajrna
    @test seqmatrix(prot, :site) == sitemajaa
    @test seqmatrix(UInt8, dna, :site) == sitemajnucint
    @test seqmatrix(UInt8, rna, :site) == sitemajnucint

    @test seqmatrix(dna, :seq) == seqmajdna
    @test seqmatrix(rna, :seq) == seqmajrna
    @test seqmatrix(prot, :seq) == seqmajaa
    @test seqmatrix(UInt8, dna, :seq) == seqmajnucint
    @test seqmatrix(UInt8, rna, :seq) == seqmajnucint

    @test seqmatrix([dna"", dna"", dna""], :site) == Matrix{DNA}(undef, (3, 0))
    @test seqmatrix([dna"", dna"", dna""], :seq) == Matrix{DNA}(undef, (0, 3))
    @test seqmatrix([rna"", rna"", rna""], :site) == Matrix{RNA}(undef, (3, 0))
    @test seqmatrix([rna"", rna"", rna""], :seq) == Matrix{RNA}(undef, (0, 3))
    @test seqmatrix(UInt8, [dna"", dna"", dna""], :site) == Matrix{UInt8}(undef, (3, 0))
    @test seqmatrix(UInt8, [dna"", dna"", dna""], :seq) == Matrix{UInt8}(undef, (0, 3))
    @test seqmatrix(UInt8, [rna"", rna"", rna""], :site) == Matrix{UInt8}(undef, (3, 0))
    @test seqmatrix(UInt8, [rna"", rna"", rna""], :seq) == Matrix{UInt8}(undef, (0, 3))

    @test_throws ArgumentError seqmatrix(dnathrow, :site)
    @test_throws ArgumentError seqmatrix(rnathrow, :seq)
    @test_throws ArgumentError seqmatrix(dna, :lol)
    @test_throws MethodError seqmatrix(AminoAcid, dna, :site)
    @test_throws ArgumentError seqmatrix(DNASequence[], :site)
    @test_throws ArgumentError seqmatrix(DNASequence[], :seq)

end
