# Note: This test suite has many hard coded values in it.
# While that is normally a bad idea, we need to guarantee reproducible results
# when using the same seed.

global SEED = 0 # Do not change this
struct MyAlphabet <: Alphabet end

@testset "Random LongSequences" begin
function test_sampler(sampler, seed, elements, firstten, T)
    rng = StableRNG(seed)
    sampled1 = [rand(rng, sampler) for i in 1:1000]
    sampled2 = rand(StableRNG(seed), sampler, 1000)
    @test sampled1 == sampled2
    @test Set(sampled1) == Set(T[convert(T, i) for i in elements])
    @test eltype(sampler) == T
    @test eltype(firstten) == T
    @test rand(StableRNG(seed), sampler, 10) == firstten
end

function test_isseq(seq, alphabettype, len)
    @test seq isa LongSequence{alphabettype}
    @test length(seq) == len
end

@testset "SamplerUniform" begin
    # Cannot instantiate with empty collection
    @test_throws ArgumentError sampler = SamplerUniform{DNA}(DNA[])

    # DNA sampler with DNA array
    sampler = SamplerUniform{DNA}([DNA_A, DNA_C, DNA_G, DNA_W])
    firstten = DNA[DNA_C, DNA_W, DNA_G, DNA_C, DNA_C, DNA_C, DNA_C, DNA_A, DNA_A, DNA_A]
    test_sampler(sampler, SEED, sampler.elems, firstten, DNA)
    # Now RNA sampler from DNA array
    sampler = SamplerUniform{RNA}([DNA_A, DNA_C, DNA_G, DNA_W])
    firstten = RNA[RNA_C, RNA_W, RNA_G, RNA_C, RNA_C, RNA_C, RNA_C, RNA_A, RNA_A, RNA_A]
    test_sampler(sampler, SEED, sampler.elems, firstten, RNA)

    # Cannot make AA sampler from DNA
    @test_throws MethodError s = SamplerUniform{AminoAcid}([DNA_A, DNA_C, DNA_G, DNA_W])

    # Automatically infer eltype
    sampler1 = SamplerUniform{DNA}([DNA_A, DNA_C])
    sampler2 = SamplerUniform([DNA_A, DNA_C])
    @test typeof(sampler2) == SamplerUniform{DNA}
    @test eltype(sampler1) == eltype(sampler2)

    # Can also infer abstract eltype
    sampler3 = SamplerUniform([DNA_A, RNA_A])
    @test eltype(sampler3) == typejoin(DNA, RNA)

    sampler = SamplerUniform{RNA}([DNA_A, DNA_C, DNA_G, DNA_W])
    @test typeof(rand(sampler)) == RNA

end # SamplerUniform

@testset "SamplerWeighted" begin
    # Must have one less weight than elements
    @test_throws ArgumentError SamplerWeighted{DNA}([DNA_A], [0.5])
    @test_throws ArgumentError SamplerWeighted{DNA}([DNA_A, DNA_C], [0.5, 0.5])
    @test_throws ArgumentError SamplerWeighted{DNA}([DNA_A], [0.5, 0.5])

    # Weights cannot exceed one
    @test_throws ArgumentError SamplerWeighted{DNA}([DNA_A, DNA_C], [1.1])
    @test_throws ArgumentError SamplerWeighted{DNA}([DNA_A, DNA_C, DNA_G], [1.00001, 0.0])

    # Weights cannot be negative
    @test_throws ArgumentError SamplerWeighted{DNA}([DNA_A, DNA_C, DNA_G], [0.5, -0.001])
    @test_throws ArgumentError SamplerWeighted{DNA}([DNA_A, DNA_C, DNA_A], [1.1, -0.2])

    # Weights always sum to one after instantiation
    s = SamplerWeighted{DNA}([DNA_A, DNA_C, DNA_G, DNA_W], [0.03, 0.2, 0.7])
    @test sum(s.probs) == 1.0
    s = SamplerWeighted{DNA}([DNA_A], [])
    @test sum(s.probs) == 1.0
    s = SamplerWeighted{DNA}([DNA_A, DNA_C], [1.0])
    @test sum(s.probs) == 1.0
    s = SamplerWeighted{DNA}([DNA_A, DNA_C, DNA_G], [0.03, 0.20000001])
    @test sum(s.probs) == 1.0

    sampler = SamplerWeighted{DNA}([DNA_N, DNA_C, DNA_W, DNA_T], [0.1, 0.2, 0.3])
    firstten = DNA[DNA_C, DNA_N, DNA_T, DNA_T, DNA_T, DNA_N, DNA_W, DNA_C, DNA_W, DNA_N]
    test_sampler(sampler, 0, sampler.elems, firstten, DNA)

    sampler = SamplerWeighted{RNA}([DNA_N, DNA_C, DNA_W, DNA_T], [0.15, 0.5, 0.2])
    firstten = RNA[RNA_C, RNA_N, RNA_U, RNA_W, RNA_U, RNA_N, RNA_C, RNA_C, RNA_C, RNA_N]
    test_sampler(sampler, 0, sampler.elems, firstten, RNA)

    @test_throws MethodError s = SamplerWeighted{AminoAcid}([DNA_A, DNA_C], [0.1])

    # Casual test to see that it can use the global RNG automatically
    sampler = SamplerWeighted{RNA}([DNA_N, DNA_C, DNA_W, DNA_T], [0.15, 0.5, 0.2])
    @test typeof(rand(sampler)) == RNA

end # SamplerWeighted

@testset "randseq Sampler" begin
    # Cannot instantiate < 0-length seq
    @test_throws ArgumentError randseq(DNAAlphabet{2}(), SamplerUniform(dna"ACG"), -1)
    @test_throws ArgumentError randseq(AminoAcidAlphabet(), SamplerUniform(aa"VTW"), -1)
    @test_throws ArgumentError randseq(RNAAlphabet{4}(), SamplerWeighted(dna"ACG", ones(2)/3), -1)

    # CAN make empty sequence with mismatching alphabets
    # Or with matching alphabets
    @test randseq(DNAAlphabet{2}(), SamplerUniform(aa"AV"), 0) == LongSequence{DNAAlphabet{2}}()
    @test randseq(DNAAlphabet{2}(), SamplerUniform(rna"U"), 2) == LongSequence{DNAAlphabet{2}}(dna"TT")

    # Cannot make nonzero sequence with mismatching alphabets
    @test_throws MethodError randseq(RNAAlphabet{4}(), SamplerUniform(aa"AV", 1))
    @test_throws MethodError randseq(AminoAcidAlphabet(), SamplerUniform(dna"ACG", 10))

    # A few samplings
    sampler = SamplerUniform(aa"TVVWYAEDK")
    @test randseq(StableRNG(SEED), AminoAcidAlphabet(), sampler, 10) == aa"VEVTKTTADY"

    sampler = SamplerUniform(dna"TGAWYKN")
    @test randseq(StableRNG(SEED), DNAAlphabet{4}(), sampler, 10) == dna"GWNGGGKTTT"

    sampler = SamplerWeighted(rna"UGCMKYN", [0.1, 0.05, 0.2, 0.15, 0.15, 0.2])
    @test randseq(StableRNG(SEED), DNAAlphabet{4}(), sampler, 10) == dna"CTNYNTKCMT"

    # Casual tests to see that it can use the global RNG automatically
    sampler = SamplerUniform(aa"TVVWYAEDK")
    seq = randseq(AminoAcidAlphabet(), sampler, 25)
    test_isseq(seq, AminoAcidAlphabet, 25)

    sampler = SamplerWeighted(rna"UGCMKYN", [0.1, 0.05, 0.2, 0.15, 0.15, 0.2])
    seq = randseq(RNAAlphabet{4}(), sampler, 25)
    test_isseq(seq, RNAAlphabet{4}, 25)

    sampler = SamplerUniform(dna"AGC")
    seq = randseq(DNAAlphabet{2}(), sampler, 25)
    test_isseq(seq, DNAAlphabet{2}, 25)
end # randseq Sampler

@testset "randseq" begin
    sampler = SamplerUniform(aa"ACDEFGHIKLMNPQRSTVWY")
    automatic = randseq(StableRNG(SEED), AminoAcidAlphabet(), 1000)
    manual = randseq(StableRNG(SEED), AminoAcidAlphabet(), sampler, 1000)
    @test automatic == manual
    @test Set(automatic) == Set(aa"ACDEFGHIKLMNPQRSTVWY")

    @test randseq(StableRNG(SEED), DNAAlphabet{4}(), 20) == dna"CTCTTCGTATGCCGTACCGT"
    @test randseq(StableRNG(SEED), DNAAlphabet{2}(), 20) == dna"CATCCGCTCCTAGGCCATTT"

    @test randseq(StableRNG(SEED), RNAAlphabet{4}(), 20) == rna"CUCUUCGUAUGCCGUACCGU"
    @test randseq(StableRNG(SEED), RNAAlphabet{2}(), 20) == rna"CAUCCGCUCCUAGGCCAUUU"

    # Casual tests to see that it can use the global RNG automatically
    seq = randseq(DNAAlphabet{2}(), 100)
    test_isseq(seq, DNAAlphabet{2}, 100)

    seq = randseq(RNAAlphabet{4}(), 100)
    test_isseq(seq, RNAAlphabet{4}, 100)

    seq = randseq(AminoAcidAlphabet(), 100)
    test_isseq(seq, AminoAcidAlphabet, 100)
end # randseq

@testset "Simple constructors" begin
    manual = randseq(StableRNG(SEED), AminoAcidAlphabet(), 100)
    automatic = randaaseq(StableRNG(SEED), 100)
    @test automatic == manual

    manual = randseq(StableRNG(SEED), DNAAlphabet{4}(), 100)
    automatic = randdnaseq(StableRNG(SEED), 100)
    @test automatic == manual

    manual = randseq(StableRNG(SEED), RNAAlphabet{4}(), 100)
    automatic = randrnaseq(StableRNG(SEED), 100)
    @test automatic == manual

    # Casual tests to see that it can use the global RNG automatically
    seq = randdnaseq(10)
    test_isseq(seq, DNAAlphabet{4}, 10)

    seq = randrnaseq(10)
    test_isseq(seq, RNAAlphabet{4}, 10)

    seq = randaaseq(10)
    test_isseq(seq, AminoAcidAlphabet, 10)
end # simple constructors

@testset "Custom alphabet" begin
    Base.length(A::MyAlphabet) = 6
    BioSequences.symbols(A::MyAlphabet) = (DNA_A, DNA_C, DNA_G, DNA_T, RNA_U, DNA_N)
    BioSequences.BitsPerSymbol(A::MyAlphabet) = BioSequences.BitsPerSymbol{8}()
    BioSequences.encode(A::MyAlphabet, x::DNA) = reinterpret(UInt8, x)
    BioSequences.encode(A::MyAlphabet, x::RNA) = reinterpret(UInt8, x) | 0x10
    Base.eltype(A::MyAlphabet) = NucleicAcid

    function BioSequences.decode(A::MyAlphabet, x::UInt64)
        uint = UInt8(x)
        if uint & 0x10 == 0x10
            return reinterpret(RNA, uint & 0x0f)
        else
            return reinterpret(DNA, uint)
        end
    end

    seq = randseq(MyAlphabet(), 1000)
    @test typeof(seq) == LongSequence{MyAlphabet}
    @test length(seq) == 1000
    @test Set(seq) == Set(symbols(MyAlphabet()))
end # Custom Alphabet
end # Entire Random LongSequences testset
