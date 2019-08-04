# NOTE: Most tests related to biological symbols are located in BioSymbols.jl.
@testset "Symbols" begin
    @testset "DNA" begin
        @test DNA_A === BioSymbols.DNA_A
        @test ACGT  === BioSymbols.ACGT
        @test ACGTN === BioSymbols.ACGTN
        @test typeof(DNA_A) === BioSymbols.DNA
    end

    @testset "RNA" begin
        @test RNA_A === BioSymbols.RNA_A
        @test ACGU  === BioSymbols.ACGU
        @test ACGUN === BioSymbols.ACGUN
        @test typeof(RNA_A) === BioSymbols.RNA
    end

    @testset "AminoAcid" begin
        @test AA_A === BioSymbols.AA_A
        @test typeof(AA_A) === BioSymbols.AminoAcid
    end

    @testset "Predicate functions" begin
        @test iscompatible(DNA_A, DNA_N)
        @test isambiguous(DNA_N)
        @test iscertain(DNA_A)
        @test isgap(DNA_Gap)
        @test ispurine(DNA_A)
        @test ispyrimidine(DNA_C)
        @test isGC(DNA_G)
    end

    @testset "Misc. functions" begin
        @test length(BioSymbols.alphabet(DNA)) == 16
        @test BioSymbols.gap(DNA) === DNA_Gap
        @test BioSymbols.complement(DNA_A) === DNA_T
    end

    @testset "Encoder" begin
        encode = BioSequences.encode
        EncodeError = BioSequences.EncodeError

        @testset "DNA" begin
            # 2 bits
            @test encode(DNAAlphabet{2}(), DNA_A) === 0x00
            @test encode(DNAAlphabet{2}(), DNA_C) === 0x01
            @test encode(DNAAlphabet{2}(), DNA_G) === 0x02
            @test encode(DNAAlphabet{2}(), DNA_T) === 0x03
            @test_throws EncodeError encode(DNAAlphabet{2}(), DNA_M)
            @test_throws EncodeError encode(DNAAlphabet{2}(), DNA_N)
            @test_throws EncodeError encode(DNAAlphabet{2}(), DNA_Gap)

            # 4 bits
            for nt in BioSymbols.alphabet(DNA)
                @test encode(DNAAlphabet{4}(), nt) === reinterpret(UInt8, nt)
            end
            @test_throws EncodeError encode(DNAAlphabet{4}(), reinterpret(DNA, 0b10000))
        end

        @testset "RNA" begin
            # 2 bits
            @test encode(RNAAlphabet{2}(), RNA_A) === 0x00
            @test encode(RNAAlphabet{2}(), RNA_C) === 0x01
            @test encode(RNAAlphabet{2}(), RNA_G) === 0x02
            @test encode(RNAAlphabet{2}(), RNA_U) === 0x03
            @test_throws EncodeError encode(RNAAlphabet{2}(), RNA_M)
            @test_throws EncodeError encode(RNAAlphabet{2}(), RNA_N)
            @test_throws EncodeError encode(RNAAlphabet{2}(), RNA_Gap)

            # 4 bits
            for nt in BioSymbols.alphabet(RNA)
                @test encode(RNAAlphabet{4}(), nt) === reinterpret(UInt8, nt)
            end
            @test_throws EncodeError encode(RNAAlphabet{4}(), reinterpret(RNA, 0b10000))
        end

        @testset "AminoAcid" begin
            @test encode(AminoAcidAlphabet(), AA_A) === 0x00
            for aa in BioSymbols.alphabet(AminoAcid)
                @test encode(AminoAcidAlphabet(), aa) === convert(UInt8, aa)
            end
            @test_throws BioSequences.EncodeError encode(AminoAcidAlphabet(), BioSymbols.AA_INVALID)
        end
    end

    @testset "Decoder" begin
        decode = BioSequences.decode
        DecodeError = BioSequences.DecodeError

        @testset "DNA" begin
            # 2 bits
            @test decode(DNAAlphabet{2}(), 0x00) === DNA_A
            @test decode(DNAAlphabet{2}(), 0x01) === DNA_C
            @test decode(DNAAlphabet{2}(), 0x02) === DNA_G
            @test decode(DNAAlphabet{2}(), 0x03) === DNA_T
            @test_throws DecodeError decode(DNAAlphabet{2}(), 0x04)
            @test_throws DecodeError decode(DNAAlphabet{2}(), 0x0e)

            # 4 bits
            for x in 0b0000:0b1111
                @test decode(DNAAlphabet{4}(), x) === reinterpret(DNA, x)
            end
            @test_throws DecodeError decode(DNAAlphabet{4}(), 0b10000)
        end

        @testset "RNA" begin
            # 2 bits
            @test decode(RNAAlphabet{2}(), 0x00) === RNA_A
            @test decode(RNAAlphabet{2}(), 0x01) === RNA_C
            @test decode(RNAAlphabet{2}(), 0x02) === RNA_G
            @test decode(RNAAlphabet{2}(), 0x03) === RNA_U
            @test_throws DecodeError decode(RNAAlphabet{2}(), 0x04)
            @test_throws DecodeError decode(RNAAlphabet{2}(), 0x0e)

            # 4 bits
            for x in 0b0000:0b1111
                @test decode(RNAAlphabet{4}(), x) === reinterpret(RNA, x)
            end
            @test_throws DecodeError decode(RNAAlphabet{4}(), 0b10000)
        end

        @testset "AminoAcid" begin
            @test decode(AminoAcidAlphabet(), 0x00) === AA_A
            for x in 0x00:0x1b
                @test decode(AminoAcidAlphabet(), x) === convert(AminoAcid, x)
            end
            @test_throws BioSequences.DecodeError decode(AminoAcidAlphabet(), 0x1c)
        end
    end
end
