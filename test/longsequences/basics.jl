@testset "Basics" begin
    @test BioSequences.has_interface(BioSequence, LongDNA{2}, [DNA_G], true)
    @test BioSequences.has_interface(BioSequence, LongDNA{4}, [DNA_G], true)
    @test BioSequences.has_interface(BioSequence, LongRNA{2}, [RNA_G], true)
    @test BioSequences.has_interface(BioSequence, LongRNA{4}, [RNA_G], true)
    @test BioSequences.has_interface(BioSequence, LongAA, [AA_G], true)

	seq = LongSequence{DNAAlphabet{2}}()
	@test isempty(seq)

	@test similar(seq) == seq

	seq = LongSequence{DNAAlphabet{2}}([DNA(i) for i in "TAGCA"])
	@test seq isa LongSequence
	@test seq isa LongSequence{DNAAlphabet{2}}
	sim = similar(seq)
	@test typeof(sim) == typeof(seq)
	@test length(sim) == length(seq)

    # Construct from other sequences
    seq = LongAA("AGCTVMN")
    @test LongSequence(seq) == LongAA("AGCTVMN")
    @test LongSequence(seq, 2:5) == LongAA("GCTV")
    @test LongSequence(SimpleSeq("AUCGU")) isa LongRNA{2}
    @test LongSequence(SimpleSeq("AUCGU")) == LongRNA{2}("AUCGU")
    LongDNA{4}(LongRNA{4}("AUCGUA")) == LongDNA{4}("ATCGTA")
end

@testset "Copy sequence" begin
    # Test copy from sequence to sequence
    function test_copy(A, str)
        seq = LongSequence{A}(str)
        str2 = String(copy(seq))
        @test (str == str2 == String(seq))
    end

    for len in [0, 1, 15, 33]
        test_copy(DNAAlphabet{4}, random_dna(len))
        test_copy(RNAAlphabet{2}, random_rna(len, [0.25, 0.25, 0.25, 0.25]))
        test_copy(AminoAcidAlphabet, random_aa(len))
    end
end # testset

@testset "Copy! sequence" begin
    function test_copy!(A, srctxt)
        src = LongSequence{A}(srctxt)
        dst = LongSequence{A}(undef, 0)
        for len in [max(0, length(src) - 3), length(src), length(src) + 4]
            resize!(dst, len)
            copy!(dst, src)
            @test String(dst) == String(src)
        end
    end

    test_copy!(DNAAlphabet{4}, random_dna(14))
    test_copy!(RNAAlphabet{2}, random_rna(55, [0.25, 0.25, 0.25, 0.25]))
    test_copy!(AminoAcidAlphabet, random_aa(18))

    #  Also works across nucleotide types!
    for N in (2,4)
        src = LongDNA{N}(random_dna(33, [0.25, 0.25, 0.25, 0.25]))
        dst = LongRNA{N}(random_rna(31, [0.25, 0.25, 0.25, 0.25]))
        copy!(dst, src)
        @test String(typeof(dst)(src)) == String(dst)
        resize!(dst, 16)
        copy!(src, dst)
        @test String(typeof(dst)(src)) == String(dst)
    end

    # Doesn't work for wrong types
    @test_throws Exception copy!(LongDNA{4}("TAG"), LongAA("WGM"))
    @test_throws Exception copy!(LongDNA{2}("TAG"), LongRNA{4}("UGM"))
end

@testset "Copyto! sequence" begin
    function test_copyto!1(A1, dst, A2, src)
        dst_ = LongSequence{A1}(dst)
        src_ = LongSequence{A2}(src)
        copyto!(dst_, src_)
        @test String(src_) == String(dst_[1:length(src_)])
    end

    test_copyto!1(DNAAlphabet{4}, random_dna(19), DNAAlphabet{4}, random_dna(17))
    test_copyto!1(RNAAlphabet{2}, random_rna(31, [0.25, 0.25, 0.25, 0.25]),
                  RNAAlphabet{2}, random_rna(11, [0.25, 0.25, 0.25, 0.25]))
    test_copyto!1(AminoAcidAlphabet, random_aa(61), AminoAcidAlphabet, random_aa(61))

    function test_copyto2!(A, F)
        for len in [10, 17, 51]
            start = rand(1:3)
            N = len - rand(1:4) - start
            src = LongSequence{A}(F(len))
            dst = LongSequence{A}(F(len))
            copyto!(dst, start, src, start + 1, N)
            @test String(dst[start:start+N-1]) == String(src[start+1:start+N])
        end
    end

    test_copyto2!(DNAAlphabet{2}, len -> random_dna(len, [0.25, 0.25, 0.25, 0.25]))
    test_copyto2!(RNAAlphabet{4}, random_rna)
    test_copyto2!(AminoAcidAlphabet, random_aa)

    # Test bug when copying to self
    src = LongDNA{4}("A"^16 * "C"^16 * "A"^16)
    copyto!(src, 17, src, 1, 32)
    @test String(src) == "A"^32 * "C"^16

    # Can't copy over edge
    dst = LongDNA{4}("TAGCA")
    @test_throws Exception copyto!(dst, 1, fill(0x61, 2), 2, 3)
end

@testset "Copy! data" begin
    function test_copy!(seq, src)
        @test String(src) == String(copy!(seq, src))
    end
    # Needed because conversion to String truncates vector.
    function test_copy!(seq, src::Vector)
        @test String(copy(src)) == String(copy!(copy(seq), src))
    end

    probs = [0.25, 0.25, 0.25, 0.25, 0.00]
    dna2 = LongDNA{2}(undef, 6)
    dna4 = LongDNA{4}(undef, 6)
    rna2 = LongRNA{2}(undef, 6)
    rna4 = LongRNA{4}(undef, 6)
    aa = LongAA(undef, 6)
    for dtype in [Vector{UInt8}, Vector{Char}, String, Test.GenericString]
        for len in [0, 1, 5, 16, 32, 100]
            test_copy!(dna2, dtype(random_dna(len, probs)))
            test_copy!(dna4, dtype(random_dna(len)))
            test_copy!(rna2, dtype(random_rna(len, probs)))
            test_copy!(rna4, dtype(random_rna(len)))
            test_copy!(aa, dtype(random_aa(len)))
        end
    end
end

@testset "Copyto! data" begin
    function test_twoarg_copyto!(seq, src)
        copyto!(seq, src)
        @test String(src[1:length(src)]) == String(src)
    end
    # Needed because conversion to String truncates vector.
    function test_twoarg_copyto!(seq, src::Vector)
        copyto!(seq, src)
        @test String(src[1:length(src)]) == String(copy(src))
    end

    probs = [0.25, 0.25, 0.25, 0.25, 0.00]
    dna2 = LongDNA{2}(undef, 50)
    dna4 = LongDNA{4}(undef, 50)
    rna2 = LongRNA{2}(undef, 50)
    rna4 = LongRNA{4}(undef, 50)
    aa = LongAA(undef, 50)
    for dtype in [Vector{UInt8}, Vector{Char}, String, Test.GenericString]
        for len in [0, 1, 10, 16, 32, 5]
            test_twoarg_copyto!(dna2, dtype(random_dna(len, probs)))
            test_twoarg_copyto!(dna4, dtype(random_dna(len)))
            test_twoarg_copyto!(rna2, dtype(random_rna(len, probs)))
            test_twoarg_copyto!(rna4, dtype(random_rna(len)))
            test_twoarg_copyto!(aa, dtype(random_aa(len)))
        end
    end

    # Five-arg copyto!
    function test_fivearg_copyto!(seq, src)
        for soff in (1, 3)
            for doff in (1, 5)
                for N in (0, 5, 18, 30)
                    copyto!(seq, doff, src, soff, N)
                    @test String(seq[doff:doff+N-1]) == String(src[soff:soff+N-1])
                end
            end
        end
    end

    for dtype in [Vector{UInt8}, Vector{Char}, String, Test.GenericString]
        test_fivearg_copyto!(dna2, dtype(random_dna(60, probs)))
        test_fivearg_copyto!(dna4, dtype(random_dna(60)))
        test_fivearg_copyto!(rna2, dtype(random_rna(60, probs)))
        test_fivearg_copyto!(rna4, dtype(random_rna(60)))
        test_fivearg_copyto!(aa, dtype(random_aa(60)))
    end

end

################

@testset "Concatenation" begin
    function test_concatenation(A, chunks)
        parts = UnitRange{Int}[]
        for i in 1:lastindex(chunks)
            start = rand(1:length(chunks[i]))
            stop = rand(start:length(chunks[i]))
            push!(parts, start:stop)
        end
        str = string([chunk[parts[i]] for (i, chunk) in enumerate(chunks)]...)
        seq = *([LongSequence{A}(chunk)[parts[i]] for (i, chunk) in enumerate(chunks)]...)
        @test String(seq) == uppercase(str)
    end

    for _ in [1, 2, 5, 10, 18, 32, 55, 64, 70]
        n = rand(1:10)
        chunks = [random_dna(rand(1:100)) for _ in 1:n]
        test_concatenation(DNAAlphabet{4}, chunks)

        chunks = [random_rna(rand(1:100)) for _ in 1:n]
        test_concatenation(RNAAlphabet{4}, chunks)

        chunks = [random_aa(rand(1:100)) for _ in 1:n]
        test_concatenation(AminoAcidAlphabet, chunks)

        probs = [0.25, 0.25, 0.25, 0.25, 0.00]
        chunks = [random_dna(rand(1:100), probs) for _ in 1:n]
        test_concatenation(DNAAlphabet{2}, chunks)

        chunks = [random_rna(rand(1:100), probs) for _ in 1:n]
        test_concatenation(RNAAlphabet{2}, chunks)
    end
end

@testset "Repetition" begin
    function test_repetition(A, chunk)
        start = rand(1:length(chunk))
        stop = rand(start:length(chunk))
        n = rand(1:10)
        str = chunk[start:stop] ^ n
        seq = LongSequence{A}(chunk)[start:stop] ^ n
        @test String(seq) == uppercase(str)
    end

    for _ in 1:10
        chunk = random_dna(rand(1:100))
        test_repetition(DNAAlphabet{4}, chunk)

        chunk = random_rna(rand(1:100))
        test_repetition(RNAAlphabet{4}, chunk)

        chunk = random_aa(rand(1:100))
        test_repetition(AminoAcidAlphabet, chunk)

        probs = [0.25, 0.25, 0.25, 0.25, 0.00]
        chunk = random_dna(rand(1:100), probs)
        test_repetition(DNAAlphabet{2}, chunk)

        chunk = random_rna(rand(1:100), probs)
        test_repetition(RNAAlphabet{2}, chunk)
    end
end

@testset "Join" begin
    function test_join(::Type{T}, seqs, result) where T
        @test join(T, seqs) == result
        @test join!(T(), seqs) == result
        long_seq = T(undef, 1000)
        @test join!(long_seq, seqs) == result
    end

    base_seq = dna"TGATGCTAVWMMKACGAS" # used for seqviews in testing
    test_join(LongDNA{4}, [dna"TACG", dna"ACCTGT", @view(base_seq[16:18])], dna"TACGACCTGTGAS")
    test_join(LongRNA{4}, Set([]), LongRNA{4}())

    base_aa_seq = aa"KMAEEHPAIYWLMN"
    test_join(LongAA, (aa"KMVLE", aa"", (@view base_aa_seq[3:6])), aa"KMVLEAEEH")
end

@testset "Length" begin
    for len in [0, 1, 2, 10, 16, 32, 1000]
        seq = LongDNA{4}(random_dna(len))
        @test length(seq) === lastindex(seq) === len

        seq = LongRNA{4}(random_rna(len))
        @test length(seq) === lastindex(seq) === len

        seq = LongAA(random_aa(len))
        @test length(seq) === lastindex(seq) === len
    end
end

@testset "Access" begin
    dna_seq = dna"ACTG"

    @test dna_seq[1] === DNA_A
    @test dna_seq[2] === DNA_C
    @test dna_seq[3] === DNA_T
    @test dna_seq[4] === DNA_G

    # Access indexes out of bounds
    @test_throws BoundsError dna_seq[-1]
    @test_throws BoundsError dna_seq[0]
    @test_throws BoundsError dna_seq[5]

    @test dna"ACTGNACTGN"[1:5] == dna"ACTGN"
    @test dna"ACTGNACTGN"[5:1] == dna""

    rna_seq = rna"ACUG"
    @test rna_seq[1] === RNA_A
    @test rna_seq[2] === RNA_C
    @test rna_seq[3] === RNA_U
    @test rna_seq[4] === RNA_G

    # Access indexes out of bounds
    @test_throws BoundsError rna_seq[-1]
    @test_throws BoundsError rna_seq[0]
    @test_throws BoundsError rna_seq[5]

    @test rna"ACUGNACUGN"[1:5] == rna"ACUGN"
    @test rna"ACUGNACUGN"[5:1] == rna""

    @test aa"KSAAV"[3] == AA_A
end

@testset "Equality" begin
    @testset "DNA" begin
        a = dna"ACTGN"
        b = dna"ACTGN"
        @test a == b
        @test dna"ACTGN" == dna"ACTGN"
        @test dna"ACTGN" != dna"ACTGA"
        @test dna"ACTGN" != dna"ACTG"
        @test dna"ACTG"  != dna"ACTGN"

        a = dna"ACGTNACGTN"
        b = dna"""
        ACGTN
        ACGTN
        """
        @test a == b
    end

    @testset "RNA" begin
        a = rna"ACUGN"
        b = rna"ACUGN"
        @test a == b
        @test rna"ACUGN" == rna"ACUGN"
        @test rna"ACUGN" != rna"ACUGA"
        @test rna"ACUGN" != rna"ACUG"
        @test rna"ACUG"  != rna"ACUGN"

        a = rna"ACUGNACUGN"
        b = rna"""
        ACUGN
        ACUGN
        """
        @test a == b
    end

    @testset "AminoAcid" begin
        a = aa"ARNDCQEGHILKMFPSTWYVX"
        b = aa"ARNDCQEGHILKMFPSTWYVX"
        @test a == b
        @test a != aa"ARNDCQEGHILKMFPSTWYXV"
        @test a != aa"ARNDCQEGHLKMFPSTWYVX"

        b = aa"""
        ARNDCQEGHI
        LKMFPSTWYV
        X
        """
        @test a == b
    end
end

@testset "Custom ASCII alphabet" begin
    @test string(LongSequence{ReducedAAAlphabet}("FSPMKH")) == "FSPMKH"
    buf = IOBuffer()
    print(buf, LongSequence{ReducedAAAlphabet}("AGTDNWLE"))
    @test String(take!(buf)) == "AGTDNWLE"

    #### Now add AsciiAlphabet capacity
    BioSequences.codetype(::ReducedAAAlphabet) = BioSequences.AsciiAlphabet()
    function BioSequences.ascii_encode(::ReducedAAAlphabet, x::UInt8)
        for sym in symbols(ReducedAAAlphabet())
            if UInt8(Char(sym)) == x
                return UInt8(encode(ReducedAAAlphabet(), sym))
            end
        end
    end

    @test string(LongSequence{ReducedAAAlphabet}("FSPMKH")) == "FSPMKH"
    buf = IOBuffer()
    print(buf, LongSequence{ReducedAAAlphabet}("AGTDNWLE"))
    @test String(take!(buf)) == "AGTDNWLE"
end