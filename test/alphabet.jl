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

"""
# Alphabets of biological symbols.

`Alphabet` is the most important type trait for `BioSequence` An `Alphabet`
represents a set of biological symbols encoded by a sequence, e.g. A, C, G
and T for a 2-bit DNA Alphabet.

* Subtypes of Alphabet are singleton structs that may or may not be parameterized.
* Alphabets span over a *finite* list of biological symbols
* The alphabet controls the encoding/decoding between a decoded element type and
    an internal data representation type.
* An `Alphabet` must never encode (using `encode`) or decode (using `decode`)
invalid data.  Other methods for check-free encoding/decoding methods may be added.

Every subtype `A` of `Alphabet` must implement:
* `Base.eltype(::Type{A})::Type{E}` for some eltype `E`, which must be a BioSymbol
* `symbols(::A)::Tuple{Vararg{E}}`. This gives tuplea of all elements of `A`.
* `encode(::A, ::E)::X` encodes an element to the internal data eltype `X`
* `decode(::A, ::X)::E` decodes an `X` to an element `E`.
* Except for `eltype` which must follow Base conventions, all functions operating
on `Alphabet` should operate on instances of the alphabet, not the type.


If you want interoperation with existing subtypes of `BioSequence`,
the element type `E` must be `UInt`, and you must also implement:
* `BitsPerSymbol(::A)::BitsPerSymbol{N}`, where the `N` must be zero
or a power of two in [1, 2, 4, 8, 16, 32, [64 for 64-bit systems]].
"""
encode = BioSequences.encode
EncodeError = BioSequences.EncodeError
decode = BioSequences.decode
DecodeError = BioSequences.DecodeError

@testset "Custom Alphabet" begin
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
        x ≥ length(DEC_LUT) && throw(DecodeError(ReducedAAAlphabet(), x))
        DEC_LUT[x + UInt(1)]
    end

    #### Tests
    @test length(symbols(ReducedAAAlphabet())) == 16
    @test all(i isa AminoAcid for i in symbols(ReducedAAAlphabet()))
    @test length(Set(symbols(ReducedAAAlphabet()))) == 16
    
    for aa in [AA_P, AA_L, AA_H, AA_M]
        data = UInt(findfirst(isequal(aa), symbols(ReducedAAAlphabet())) - 1)
        @test encode(ReducedAAAlphabet(), aa) === data
        @test decode(ReducedAAAlphabet(), data) === aa
    end
    @test_throws EncodeError encode(ReducedAAAlphabet(), AA_V)
    @test_throws EncodeError encode(ReducedAAAlphabet(), AA_I)
    @test_throws EncodeError encode(ReducedAAAlphabet(), AA_R)
    @test_throws EncodeError encode(ReducedAAAlphabet(), AA_Gap)
    @test_throws EncodeError encode(ReducedAAAlphabet(), reinterpret(AminoAcid, 0xff))

    @test_throws DecodeError decode(ReducedAAAlphabet(), UInt(16))
    @test_throws DecodeError decode(ReducedAAAlphabet(), UInt(255))
    @test_throws DecodeError decode(ReducedAAAlphabet(), UInt(432881))
    @test_throws DecodeError decode(ReducedAAAlphabet(), typemax(UInt))
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

        # 4 bits
        for nt in BioSymbols.alphabet(DNA)
            @test encode(DNAAlphabet{4}(), nt) === UInt(reinterpret(UInt8, nt))
        end
        @test_throws EncodeError encode(DNAAlphabet{4}(), reinterpret(DNA, 0b10000))
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

        # 4 bits
        for nt in BioSymbols.alphabet(RNA)
            @test encode(RNAAlphabet{4}(), nt) === UInt(reinterpret(UInt8, nt))
        end
        @test_throws EncodeError encode(RNAAlphabet{4}(), reinterpret(RNA, 0b10000))
    end

    @testset "AminoAcid" begin
        @test encode(AminoAcidAlphabet(), AA_A) === UInt(0x00)
        for aa in BioSymbols.alphabet(AminoAcid)
            @test encode(AminoAcidAlphabet(), aa) === convert(UInt, aa)
        end
        @test_throws BioSequences.EncodeError encode(AminoAcidAlphabet(), BioSymbols.AA_INVALID)
    end
end

@testset "Decoding DNA/RNA/AminoAcid" begin
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
