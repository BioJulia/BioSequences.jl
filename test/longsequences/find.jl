@testset "Find" begin
    seq = dna"TACT-WKKATYACCG"

    @test findfirst(==(DNA_T), seq) == 1
    @test findfirst(==(DNA_A), seq) == 2
    @test findfirst(isequal(DNA_C), seq) == 3
    @test findfirst(==(DNA_Gap), seq) == 5
    @test findfirst(isequal(DNA_W), seq) == 6
    @test findfirst(==(DNA_Y), seq) == 11
    @test findfirst(==(DNA_G), seq) == 15
    @test findfirst(==(DNA_M), seq) === nothing
    @test findfirst(==(RNA_A), seq) === nothing
    @test findfirst(==(AA_W), seq) === nothing
    @test findfirst(==('C'), seq) === nothing

    @test findlast(isequal(DNA_C), seq) == 14
    @test findlast(isequal(DNA_G), seq) == 15
    @test findlast(isequal(DNA_Y), seq) == 11
    @test findlast(isequal(DNA_W), seq) == 6
    @test findlast(isequal(RNA_Y), seq) === nothing
    @test findlast(isequal(RNA_Gap), seq) === nothing

    #               1   5    10   15   20
    seq = view(aa"WKRRLY-JNMOPQALKAVWYOPQRT", 3:22)
    @test findnext(==(AA_R), seq, 1) == 1
    @test findnext(==(AA_R), seq, 3) === nothing
    @test findnext(==(AA_J), seq, 6) == 6
    @test findnext(==(AA_J), seq, 7) === nothing
    @test findnext(==(AA_Q), seq, 5) == 11
    @test findnext(==(AA_V), seq, 1) == 16
    @test findnext(==(AA_Y), seq, 4) == 4
    @test findnext(==(AA_Y), seq, 5) == 18
    @test findnext(==(DNA_M), seq, 1) === nothing
    @test findnext(==(RNA_A), seq, 1) === nothing
    @test findnext(==(AA_C), seq, 1) === nothing
    @test findnext(==(19), seq, 1) === nothing

    @test findprev(==(AA_R), seq, 20) == 2
    @test findprev(==(AA_K), seq, 15) == 14
    @test findprev(==(AA_Y), seq, 20) == 18
    @test findprev(==(AA_Y), seq, 17) == 4
    @test findlast(==(AA_Q), seq) == 11
    @test findprev(==(AA_J), seq, 16) == 6
    @test findprev(==(AA_Gap), seq, 9) == 5
    @test findprev(==(AA_O), seq, 3) === nothing
    @test findprev(==(AA_L), seq, 2) === nothing
    @test findprev(==(AA_H), seq, 17) === nothing

    seq = LongRNA{2}(rna"UAGUCGUAGC")
    @test findfirst(==(RNA_U), seq) == 1
    @test findfirst(==(RNA_A), seq) == 2
    @test findfirst(==(RNA_G), seq) == 3
    @test findfirst(==(RNA_C), seq) == 5
    @test findfirst(==(RNA_Gap), seq) === nothing
    @test findlast(==(RNA_G), seq) == 9

    # View with head and chunk only
    #                        1   5    10   15   20 
    seq = view(dna"TAGTCKYDDBN--ACGKMYAGCTATAYYKKM-", 11:32)
    @test findnext(==(DNA_Gap), seq, 1) == 2
    @test findnext(==(DNA_C), seq, 2) == 5
    @test findnext(==(DNA_Y), seq, 2) == 9
    @test findnext(==(DNA_Y), seq, 7) == 9
    @test findnext(==(DNA_M), seq, 8) == 8
    @test findnext(==(DNA_M), seq, 9) == 21

    # Very large sequence so it hits the inner SIMD loop
    seq = randdnaseq(541)
    seq[100] = DNA_Gap
    @test findfirst(==(DNA_Gap), seq) == findfirst(isgap, seq) == 100
    seq[100] = DNA_A
    seq[500] = DNA_Gap
    @test findlast(==(DNA_Gap), seq) == findlast(isgap, seq) == 500
    seq[500] = DNA_A
    @test findfirst(==(DNA_Gap), seq) == findlast(isgap, seq) == nothing

    # View with only tail
    #               1234
    seq = view(aa"KWYPAV-L", 3:6)
    @test findnext(==(AA_K), seq, 1) === nothing
    @test findnext(==(AA_Y), seq, 1) == 1
    @test findnext(==(AA_V), seq, 2) == 4
    @test findnext(==(AA_A), seq, 3) == 3
    @test findnext(==(AA_A), seq, 4) === nothing
    @test findnext(==(RNA_Y), seq, 1) === nothing

    @test findprev(==(AA_K), seq, 4) === nothing
    @test findprev(==(AA_P), seq, 3) == 2
    @test findprev(==(AA_A), seq, 2) === nothing
    @test findprev(==(AA_Y), seq, 4) == 1

    # Empty view
    seq = view(dna"ATCGGACTGTAGTATGAACGGATATATGCTTATAATG", 19:18)
    for i in [DNA_A, DNA_C, DNA_G, DNA_T, DNA_Gap]
        @test findfirst(==(i), seq) === nothing
        @test findlast(==(i), seq) === nothing
    end
end

@testset "SIMD find" begin
    # We exploit that all the following methods make use of the same
    # functions.
    # So, we first thoroughly check one of the methods for correctness,
    # and then we only need to check that the bitwise kernel of the rest is correct

    # The thorough one - test the right indices are found, given that the kernel
    # is correct
    @testset "isambiguous 4bit" begin
        seq = randdnaseq(400)
        indices = [1, 2, 5, 11, 21, 50, 211, 380, 391, 399, 400]
        @test isnothing(findfirst(isambiguous, seq))
        for i in indices
            seq[i] = DNA_W
        end
        for i in indices
            @test findfirst(isambiguous, seq) == i
            seq[i] = DNA_A
        end
        for i in indices
            seq[i] = DNA_K
        end
        for i in reverse(indices)
            @test findlast(isambiguous, seq) == i
            seq[i] = DNA_T
        end
        @test isnothing(findfirst(isambiguous, seq))

        # Some random tests
        @test findnext(isambiguous, dna"ATGCTGTA", 2) === nothing

        seq = dna"ATSCTGCA"
        @test findnext(isambiguous, seq, 2) === 3
        @test findnext(isambiguous, seq, 3) === 3
        @test findnext(isambiguous, seq, 4) === nothing
        @test findprev(isambiguous, seq, 8) === 3
        @test findprev(isambiguous, seq, 3) === 3
        @test findprev(isambiguous, seq, 2) === nothing
        @test_throws BoundsError findnext(isambiguous, seq, 0)
        @test_throws BoundsError findprev(isambiguous, seq, 10)

        seq = dna"ATSCTGMA"
        @test findnext(isambiguous, seq, 3) === 3
        @test findnext(isambiguous, seq, 4) === 7
        @test findnext(isambiguous, seq, 8) === nothing
        @test findprev(isambiguous, seq, 8) === 7
        @test findprev(isambiguous, seq, 7) === 7
        @test findprev(isambiguous, seq, 6) === 3
        @test findprev(isambiguous, seq, 2) === nothing
    end

    @testset "Various kernels" begin
        for (f, A) in Any[
            (isambiguous, DNAAlphabet{4}()),
            (isambiguous, AminoAcidAlphabet()),
            (!isambiguous, DNAAlphabet{4}()),
            (!isambiguous, AminoAcidAlphabet()),
            (iscertain, DNAAlphabet{4}()),
            (!iscertain, AminoAcidAlphabet()),
            (isgap, DNAAlphabet{4}()),
            (!isgap, AminoAcidAlphabet()),
        ]
            yes, no = Any[], Any[]
            for i in A
                push!(f(i) ? yes : no, i)
            end
            seq = LongSequence{typeof(A)}(undef, 3)
            for i in eachindex(seq)
                seq[i] = first(no)
            end

            @test isnothing(findfirst(f, seq))
            @test isnothing(findlast(f, seq))

            # Test validity of all 'yes'
            for i in yes
                seq[2] = i
                @test findfirst(f, seq) === 2
                @test findlast(f, seq) === 2
            end

            # Test validity of all 'no'
            for i in no
                seq[1] = i
                seq[3] = i
                @test findfirst(f, seq) === 2
                @test findlast(f, seq) === 2
            end
        end
    end

    @testset "Noop find methods" begin
        for s in Any[randseq(DNAAlphabet{2}(), 256), randseq(RNAAlphabet{2}(), 256), SimpleSeq(randdnaseq(256))]
            @test isnothing(findfirst(isgap, s))
            @test isnothing(findlast(isgap, s))

            @test isnothing(findfirst(isambiguous, s))
            @test isnothing(findlast(isambiguous, s))

            @test isnothing(findfirst(!iscertain, s))
            @test isnothing(findlast(!iscertain, s))

            for f in [iscertain, !isambiguous, !isgap]
                for i in [10, 41, 256, 1]
                    @test findnext(f, s, i) === i
                    @test findprev(f, s, i) === i
                    @test_throws BoundsError findnext(f, s, 0)
                    @test_throws BoundsError findprev(f, s, lastindex(s) + 1)
                    @test isnothing(findnext(f, s, lastindex(s) + 1))
                    @test isnothing(findprev(f, s, 0))
                end
            end
        end
    end
end

@testset "Search" begin
    #         0000000001111
    #         1234567890123
    seq = dna"NNNNNGATCGATC"

    # Check default search.
    @test findall(DNA_A, seq) == [7, 11]
    @test findall(ExactSearchQuery(dna"A"), seq) == [7:7, 11:11]

    # Check overlap key argument.
    @test findall(ExactSearchQuery(dna"GATC", iscompatible), seq; overlap = false) == [1:4, 6:9, 10:13]
    @test findall(ExactSearchQuery(dna"GATC", iscompatible), seq; overlap = true) == [1:4, 2:5, 6:9, 10:13]

    # Check mapping of indices.
    @test findall(DNA_A, seq, 7:11) == [7, 11]
    @test findall(ExactSearchQuery(dna"A"), seq, 7:11) == [7:7, 11:11]

    # Check empty return type.
    @test findall(DNA_A, dna"GGGG") |> typeof == Vector{Int}
    @test findall(ExactSearchQuery(dna"A"), dna"GGGG") |> typeof == Vector{UnitRange{Int}}

    @test findall(isequal(DNA_A), dna"ACGTAC") == [1, 5]
    @test findall(i -> true, aa"ACGTA") == collect(1:5)
    @test findall(i -> true, aa"") == Int[]
    @test findall(i -> i == AA_A, rna"AGCA") == Int[]
end
