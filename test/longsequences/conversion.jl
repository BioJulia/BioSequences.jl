@testset "Constructing empty sequences" begin
    @test isempty(LongDNA{4}())
    @test isempty(LongRNA{4}())
    @test isempty(LongAA())
end

@testset "Constructing uninitialized sequences" begin
    @test isa(LongSequence{DNAAlphabet{2}}(undef, 0), LongSequence)
    @test isa(LongSequence{DNAAlphabet{4}}(undef, 10), LongSequence)
    @test isa(LongSequence{RNAAlphabet{2}}(undef, 0), LongSequence)
    @test isa(LongSequence{RNAAlphabet{4}}(undef, 10), LongSequence)
    @test isa(LongSequence{AminoAcidAlphabet}(undef, 10), LongSequence)

    @test_throws ArgumentError LongSequence{DNAAlphabet{2}}(undef, -1)
    @test_throws ArgumentError LongSequence{DNAAlphabet{4}}(undef, -1)
    @test_throws ArgumentError LongSequence{RNAAlphabet{2}}(undef, -1)
    @test_throws ArgumentError LongSequence{RNAAlphabet{4}}(undef, -1)
    @test_throws ArgumentError LongSequence{AminoAcidAlphabet}(undef, -1)
end

@testset "Conversion from/to strings" begin
    # Check that sequences in strings survive round trip conversion:
    #   String → LongSequence → String
    function test_string_construction(A::Type, seq::AbstractString)
        @test String(LongSequence{A}(seq)) == uppercase(seq)
    end

    function test_string_parse(A::Type, seq::AbstractString)
        @test parse(LongSequence{A}, seq) == LongSequence{A}(seq)
    end

    for len in [0, 1, 3, 10, 32, 100]
        test_string_construction(DNAAlphabet{4}, random_dna(len))
        test_string_construction(DNAAlphabet{4}, SubString(random_dna(len), 1:len))
        test_string_construction(DNAAlphabet{4}, lowercase(random_dna(len)))
        test_string_construction(RNAAlphabet{4}, lowercase(random_rna(len)))
        test_string_construction(RNAAlphabet{4}, random_rna(len))
        test_string_construction(AminoAcidAlphabet, random_aa(len))
        test_string_construction(AminoAcidAlphabet, lowercase(random_aa(len)))

        test_string_parse(DNAAlphabet{4}, random_dna(len))
        test_string_parse(DNAAlphabet{4}, SubString(random_dna(len), 1:len))
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
    @test isa(dna"ACGTMRWSYKVHDBN-", LongDNA{4})
    @test isa(rna"ACGUMRWSYKVHDBN-", LongRNA{4})
    @test isa(aa"ARNDCQEGHILKMFPSTWYVBJZXOU*-", LongAA)

    # Non-nucleotide characters should throw
    @test_throws Exception LongDNA{4}("ACCNNCATTTTTTAGATXATAG")
    @test_throws Exception LongRNA{4}("ACCNNCATTTTTTAGATXATAG")
    @test_throws Exception LongAA("ATGHLMY@ZACAGNM")
end

@testset "Construction from vectors" begin
    function test_vector_construction(A, seq::AbstractString)
        T = eltype(A)
        xs = T[convert(T, c) for c in seq]
        @test LongSequence{A}(xs) == LongSequence{A}(seq)
    end

    # Construct from abstract vector
    LongDNA{4}(0x61:0x64) == LongDNA{4}("ABCD")
    LongDNA{4}(0x61:0x64, 3:4) == LongDNA{4}("CD")
    LongRNA{2}(0x61:0x61) == LongRNA{2}("A")
    LongDNA{4}(Test.GenericString("AGCTMYWK")) == LongDNA{4}("AGCTMYWK")
    LongAA(Test.GenericString("KMSPIYT")) == LongAA("KMSPIYT")

    for len in [0, 1, 10, 32, 1000]
        test_vector_construction(DNAAlphabet{4}, random_dna(len))
        test_vector_construction(RNAAlphabet{4}, random_rna(len))
        test_vector_construction(AminoAcidAlphabet, random_aa(len))

        probs = [0.25, 0.25, 0.25, 0.25, 0.00]
        test_vector_construction(DNAAlphabet{2}, random_dna(len, probs))
        test_vector_construction(RNAAlphabet{2}, random_rna(len, probs))
    end
end

@testset "Encode_copy!" begin
    # Note: Other packages use this function, so we need to test it
    # Even though this is NOT exported or part of the API in a normal sense
    function test_copyto!(dst::LongSequence, doff, src, soff, N)
        BioSequences.copyto!(dst, doff, src, soff, N)
        @test String(dst[doff:doff+N-1]) == String(src[soff:soff+N-1])
    end

    probs = [0.25, 0.25, 0.25, 0.25]
    for len in [0, 1, 10, 100]
        for f in [identity, Vector{Char}, Vector{UInt8}]
            test_copyto!(LongSequence{DNAAlphabet{2}}(undef, len), 1, f(random_dna(len, probs)), 1, len)
            test_copyto!(LongSequence{RNAAlphabet{2}}(undef, len), 1, f(random_rna(len, probs)), 1, len)
            test_copyto!(LongSequence{DNAAlphabet{4}}(undef, len), 1, f(random_dna(len)), 1, len)
            test_copyto!(LongSequence{RNAAlphabet{4}}(undef, len), 1, f(random_rna(len)), 1, len)
            test_copyto!(LongSequence{AminoAcidAlphabet}(undef, len), 1, f(random_aa(len)), 1, len)
        end
    end

    for len in [10, 32, 100]
        for f in [identity, Vector{Char}, Vector{UInt8}]
            test_copyto!(LongSequence{DNAAlphabet{2}}(undef, len+7), 5, f(random_dna(len+11, probs)), 3, len)
            test_copyto!(LongSequence{RNAAlphabet{2}}(undef, len+7), 5, f(random_rna(len+11, probs)), 3, len)
            test_copyto!(LongSequence{DNAAlphabet{4}}(undef, len+7), 5, f(random_dna(len+11)), 3, len)
            test_copyto!(LongSequence{RNAAlphabet{4}}(undef, len+7), 5, f(random_rna(len+11)), 3, len)
            test_copyto!(LongSequence{AminoAcidAlphabet}(undef, len+7), 5, f(random_aa(len+11)), 3, len)
        end
    end
end

@testset "Convert to same type" begin
	function test_same_conversion(seq)
		@test convert(typeof(seq), seq) === seq
	end

	test_same_conversion(random_dna(20))
	test_same_conversion(random_rna(20))
	test_same_conversion(random_aa(20))
end

@testset "Conversion between 2-bit and 4-bit encodings" begin
    function test_conversion(A1, A2, seq)
        @test convert(LongSequence{A1}, LongSequence{A2}(seq)) == LongSequence{A1}(seq)
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
    @test_throws EncodeError convert(LongSequence{DNAAlphabet{2}}, dna"AN")
    @test_throws EncodeError convert(LongSequence{RNAAlphabet{2}}, rna"AN")

    # test promotion
    a = LongSequence{DNAAlphabet{2}}("ATCG")
    b = LongSequence{DNAAlphabet{4}}("ATCG")
    c = LongSequence{RNAAlphabet{2}}("AUCG")
    d = LongSequence{RNAAlphabet{4}}("AUCG")

    @test typeof(promote(a, b)) == Tuple{LongSequence{DNAAlphabet{4}},LongSequence{DNAAlphabet{4}}}
    @test typeof(promote(c, d)) == Tuple{LongSequence{RNAAlphabet{4}},LongSequence{RNAAlphabet{4}}}
    @test_throws ErrorException typeof(promote(a, d))
    @test_throws ErrorException typeof(promote(a, b, d))
end

@testset "Conversion between RNA and DNA" begin
    @test convert(LongRNA{4}, LongDNA{4}("ACGTN")) == rna"ACGUN"
    @test convert(LongDNA{4}, LongRNA{4}("ACGUN")) == dna"ACGTN"
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
    @test_throws ArgumentError seqmatrix(LongDNA{4}[], :site)
    @test_throws ArgumentError seqmatrix(LongDNA{4}[], :seq)

end
