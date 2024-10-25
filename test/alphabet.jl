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
end

@testset "Basic" begin
    @test length(DNAAlphabet{4}()) == 16
    @test length(RNAAlphabet{2}()) == 4
    @test length(AminoAcidAlphabet()) == 28

    @test BioSequences.symbols(RNAAlphabet{2}()) == (RNA_A, RNA_C, RNA_G, RNA_U)
end

encode = BioSequences.encode
tryencode = BioSequences.tryencode
EncodeError = BioSequences.EncodeError
decode = BioSequences.decode

@testset "EncodeError" begin
    @test_throws EncodeError encode(DNAAlphabet{4}(), RNA_U)
    @test_throws EncodeError encode(DNAAlphabet{2}(), DNA_M)
    @test_throws EncodeError encode(DNAAlphabet{4}(), AA_C)
    @test_throws EncodeError encode(AminoAcidAlphabet(), DNA_C)
    @test_throws EncodeError encode(AminoAcidAlphabet(), RNA_N)
    @test_throws EncodeError encode(RNAAlphabet{2}(), DNA_C)
    @test_throws EncodeError encode(RNAAlphabet{2}(), RNA_K)
end

# NOTE: See the docs for the interface of Alphabet
struct ReducedAAAlphabet <: Alphabet end

Base.eltype(::Type{ReducedAAAlphabet}) = AminoAcid
BioSequences.BitsPerSymbol(::ReducedAAAlphabet) = BioSequences.BitsPerSymbol{4}()
function BioSequences.symbols(::ReducedAAAlphabet)
    (AA_L, AA_C, AA_A, AA_G, AA_S, AA_T, AA_P, AA_F,
    AA_W, AA_E, AA_D, AA_N, AA_Q, AA_K, AA_H, AA_M)
end

(ENC_LUT, DEC_LUT) = let
    enc_lut = fill(0xff, length(alphabet(AminoAcid)))
    dec_lut = fill(AA_A, length(symbols(ReducedAAAlphabet())))
    for (i, aa) in enumerate(symbols(ReducedAAAlphabet()))
        enc_lut[reinterpret(UInt8, aa) + 0x01] = i - 1
        dec_lut[i] = aa
    end
    (Tuple(enc_lut), Tuple(dec_lut))
end

function BioSequences.encode(::ReducedAAAlphabet, aa::AminoAcid)
    i = reinterpret(UInt8, aa)
    (i ≥ length(ENC_LUT) || ENC_LUT[i + 0x01] === 0xff) && throw(EncodeError(ReducedAAAlphabet(), aa))
    ENC_LUT[i + 0x01] % UInt
end

function BioSequences.decode(::ReducedAAAlphabet, x::UInt)
    DEC_LUT[x + UInt(1)]
end

@testset "Custom Alphabet" begin
    @test BioSequences.has_interface(Alphabet, ReducedAAAlphabet())
    @test length(symbols(ReducedAAAlphabet())) == 16
    @test all(i isa AminoAcid for i in symbols(ReducedAAAlphabet()))
    @test length(Set(symbols(ReducedAAAlphabet()))) == 16
    
    for aa in [AA_P, AA_L, AA_H, AA_M]
        data = UInt(findfirst(isequal(aa), symbols(ReducedAAAlphabet())) - 1)
        @test encode(ReducedAAAlphabet(), aa) === data
        @test decode(ReducedAAAlphabet(), data) === aa
    end

    str = "NSTPHML"
    @test String(LongSequence{ReducedAAAlphabet}(str)) == str

    @test_throws EncodeError encode(ReducedAAAlphabet(), AA_V)
    @test_throws EncodeError encode(ReducedAAAlphabet(), AA_I)
    @test_throws EncodeError encode(ReducedAAAlphabet(), AA_R)
    @test_throws EncodeError encode(ReducedAAAlphabet(), AA_Gap)
    @test_throws EncodeError encode(ReducedAAAlphabet(), reinterpret(AminoAcid, 0xff))
end

@testset "Encoding DNA/RNA/AminoAcid" begin
    @testset "DNA" begin
        # 2 bits
        @test encode(DNAAlphabet{2}(), DNA_A) === UInt(0x00)
        @test encode(DNAAlphabet{2}(), DNA_C) === UInt(0x01)
        @test encode(DNAAlphabet{2}(), DNA_G) === UInt(0x02)
        @test encode(DNAAlphabet{2}(), DNA_T) === UInt(0x03)
        @test_throws EncodeError encode(DNAAlphabet{2}(), DNA_M)
        @test_throws EncodeError encode(DNAAlphabet{2}(), DNA_N)
        @test_throws EncodeError encode(DNAAlphabet{2}(), DNA_Gap)

        @test tryencode(DNAAlphabet{2}(), DNA_A) == UInt(0x00)
        @test tryencode(DNAAlphabet{2}(), DNA_C) == UInt(0x01)
        @test tryencode(DNAAlphabet{2}(), DNA_G) == UInt(0x02)
        @test tryencode(DNAAlphabet{2}(), DNA_T) == UInt(0x03)
        @test tryencode(DNAAlphabet{2}(), DNA_M) === nothing
        @test tryencode(DNAAlphabet{2}(), DNA_N) === nothing
        @test tryencode(DNAAlphabet{2}(), DNA_Gap) === nothing
        @test tryencode(DNAAlphabet{2}(), RNA_G) === nothing

        # 4 bits
        for nt in BioSymbols.alphabet(DNA)
            @test encode(DNAAlphabet{4}(), nt) === UInt(reinterpret(UInt8, nt))
        end
    end

    @testset "RNA" begin
        # 2 bits
        @test encode(RNAAlphabet{2}(), RNA_A) === UInt(0x00)
        @test encode(RNAAlphabet{2}(), RNA_C) === UInt(0x01)
        @test encode(RNAAlphabet{2}(), RNA_G) === UInt(0x02)
        @test encode(RNAAlphabet{2}(), RNA_U) === UInt(0x03)
        @test_throws EncodeError encode(RNAAlphabet{2}(), RNA_M)
        @test_throws EncodeError encode(RNAAlphabet{2}(), RNA_N)
        @test_throws EncodeError encode(RNAAlphabet{2}(), RNA_Gap)

        @test tryencode(RNAAlphabet{2}(), RNA_A) == UInt(0x00)
        @test tryencode(RNAAlphabet{2}(), RNA_C) == UInt(0x01)
        @test tryencode(RNAAlphabet{2}(), RNA_G) == UInt(0x02)
        @test tryencode(RNAAlphabet{2}(), RNA_U) == UInt(0x03)
        @test tryencode(RNAAlphabet{2}(), RNA_M) === nothing
        @test tryencode(RNAAlphabet{2}(), RNA_N) === nothing
        @test tryencode(RNAAlphabet{2}(), RNA_Gap) === nothing
        @test tryencode(RNAAlphabet{2}(), DNA_G) === nothing

        # 4 bits
        for nt in BioSymbols.alphabet(RNA)
            @test encode(RNAAlphabet{4}(), nt) === UInt(reinterpret(UInt8, nt))
            @test tryencode(RNAAlphabet{4}(), nt) ===  UInt(reinterpret(UInt8, nt))
        end
    end

    @testset "AminoAcid" begin
        @test encode(AminoAcidAlphabet(), AA_A) === UInt(0x00)
        for aa in BioSymbols.alphabet(AminoAcid)
            @test encode(AminoAcidAlphabet(), aa) === convert(UInt, reinterpret(UInt8, aa))
            @test tryencode(AminoAcidAlphabet(), aa) === convert(UInt, reinterpret(UInt8, aa))
        end
        @test_throws BioSequences.EncodeError encode(AminoAcidAlphabet(), BioSymbols.AA_INVALID)
        @test tryencode(AminoAcidAlphabet(), reinterpret(AminoAcid, typemax(UInt8))) === nothing
        @test tryencode(AminoAcidAlphabet(), BioSymbols.AA_INVALID) === nothing
    end
end

@testset "Decoding DNA/RNA/AminoAcid" begin
    @testset "DNA" begin
        # 2 bits
        @test decode(DNAAlphabet{2}(), 0x00) === DNA_A
        @test decode(DNAAlphabet{2}(), 0x01) === DNA_C
        @test decode(DNAAlphabet{2}(), 0x02) === DNA_G
        @test decode(DNAAlphabet{2}(), 0x03) === DNA_T

        # 4 bits
        for x in 0b0000:0b1111
            @test decode(DNAAlphabet{4}(), x) === reinterpret(DNA, x)
        end
    end

    @testset "RNA" begin
        # 2 bits
        @test decode(RNAAlphabet{2}(), 0x00) === RNA_A
        @test decode(RNAAlphabet{2}(), 0x01) === RNA_C
        @test decode(RNAAlphabet{2}(), 0x02) === RNA_G
        @test decode(RNAAlphabet{2}(), 0x03) === RNA_U

        # 4 bits
        for x in 0b0000:0b1111
            @test decode(RNAAlphabet{4}(), x) === reinterpret(RNA, x)
        end
    end

    @testset "AminoAcid" begin
        @test decode(AminoAcidAlphabet(), 0x00) === AA_A
        for x in 0x00:0x1b
            @test decode(AminoAcidAlphabet(), x) === reinterpret(AminoAcid, x)
        end
    end
end

@testset "Interface" begin
    @test BioSequences.has_interface(Alphabet, DNAAlphabet{2}())
    @test BioSequences.has_interface(Alphabet, DNAAlphabet{4}())
    @test BioSequences.has_interface(Alphabet, RNAAlphabet{2}())
    @test BioSequences.has_interface(Alphabet, RNAAlphabet{4}())
    @test BioSequences.has_interface(Alphabet, AminoAcidAlphabet())
end

@testset "Guess parsing and guess alphabet" begin
    for (A, Ss) in [
        (DNAAlphabet{2}(), ["", "TAGCT", "AAA"]),
        (RNAAlphabet{2}(), ["U", "CGGUAG", "CCCCU"]),
        (DNAAlphabet{4}(), ["W", "HKW--", "TAGCTATSG", "TAGC-TAG"]),
        (RNAAlphabet{4}(), ["WU", "HAUH-CD", "VKMNSU"]),
        (AminoAcidAlphabet(), ["Q", "PLBMW", "***"])
    ]
        for S in Ss
            for T in [String, SubString, Vector{UInt8}, Test.GenericString]
                @test guess_alphabet(T(S)) == A
                @test bioseq(T(S)) isa LongSequence{typeof(A)}
            end
        end
    end
    for (index, S) in [
        (4, "QMN!KK"),
        (7, "ATCGAT???"),
        (1, ";")
    ]
        for T in [String, SubString, Vector{UInt8}, Test.GenericString]
            @test guess_alphabet(T(S)) == index
            @test_throws BioSequences.EncodeError bioseq(T(S))
        end
    end
end
