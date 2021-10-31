
@testset "Getindex" begin
    @testset "Scalar getindex" begin
        seq = SimpleSeq([RNA_C, RNA_G, RNA_C, RNA_A])
        @test_throws BoundsError seq[0]
        @test_throws BoundsError seq[5]

        @test first(seq) == RNA_C
        @test last(seq) == RNA_A
        
        @test seq[2] == RNA_G
        @test seq[3] == RNA_C
    end

    @testset "Getindex w. bool array" begin
        char_arr = map(RNA, collect("ACGUAGCUAUUAUACCC"))
        seq = SimpleSeq(char_arr) # len = 17
        @test_throws BoundsError seq[fill(true, length(char_arr) - 1)]
        @test_throws BoundsError seq[fill(true, length(char_arr) + 1)]

        for i in 1:3
            arr = rand(Bool, length(seq))
            seq2 = seq[arr]
            @test length(seq2) == count(arr)
            @test seq2 == SimpleSeq(char_arr[arr])
        end
    end

    @testset "Getindex w. unit range" begin
        char_arr = map(RNA, collect("ACGUAGAUUAUUUCUCCAA")) # len = 19
        seq = SimpleSeq(char_arr)

        @test seq[1:0] == empty(seq)
        @test seq[5:4] == empty(seq)

        @test_throws BoundsError seq[-1:2]
        @test_throws BoundsError seq[5:22]
        @test_throws BoundsError seq[0:19]

        @test seq[2:5] == SimpleSeq(map(RNA, collect("CGUA")))
        @test seq[1:9] == SimpleSeq(map(RNA, collect("ACGUAGAUU")))
        @test seq[1:18] == SimpleSeq(map(RNA, collect("ACGUAGAUUAUUUCUCCA")))
        @test seq[1:19] == seq
    end

    @testset "Getindex w. colon" begin
        char_arr = map(RNA, collect("ACGUAGAUUAUUUCUCCAA")) # len = 19
        seq = SimpleSeq(char_arr)
        seq[:] == seq
        seq[:] !== seq

        seq = SimpleSeq([RNA_A, RNA_C])
        seq[:] == SimpleSeq([RNA_A, RNA_C])
        seq[:] !== seq

        seq = empty(seq)
        seq[:] == empty(seq)
        seq[:] !== empty(seq)
    end

    @testset "Getindex w. integer array" begin
        char_arr = map(RNA, collect("AGCGUAUAGCGA")) # len = 12
        seq = SimpleSeq(char_arr)

        seq[Int[]] == empty(seq)
        @test_throws BoundsError seq[[-1, 2, 4]]
        @test_throws BoundsError seq[[3, 7, 9, 0]]
        @test_throws BoundsError seq[[5, 3, 13]]

        seq[[5, 2, 1]] == SimpleSeq([RNA_U, RNA_C, RNA_A])
        seq[1:2:11] == SimpleSeq(seq.x[1:2:11])
    end
end

@testset "Setindex!" begin
    @testset "Scalar setindex!" begin
        seq = SimpleSeq([RNA_U, RNA_U, RNA_A])

        @test_throws BoundsError seq[0] = RNA_A
        @test_throws BoundsError seq[-3] = RNA_A
        @test_throws BoundsError seq[4] = RNA_U

        @test_throws MethodError seq[1] = AA_A
        @test_throws MethodError seq[1] = UInt8(1)

        seq[2] = RNA_U
        @test seq == SimpleSeq([RNA_U, RNA_U, RNA_A])

        seq[1] = RNA_G
        @test seq == SimpleSeq([RNA_G, RNA_U, RNA_A])

        # RNA/DNA can be converted freely
        seq[3] = DNA_T
        @test seq == SimpleSeq([RNA_G, RNA_U, RNA_U])
    end

    @testset "Setindex w. bool array" begin
        function test_bool_arr(s::SimpleSeq, mask::AbstractArray{Bool})
            n = count(mask)
            seq2 = random_simple(n)
            cp = copy(s)
            @test s[mask] == SimpleSeq(cp.x[mask])
        end

        random_simple(len::Integer) = SimpleSeq(rand([RNA_A, RNA_C, RNA_G, RNA_U], len))
        seq = random_simple(19)

        @test_throws DimensionMismatch seq[trues(19)] = random_simple(18)
        @test_throws BoundsError seq[trues(18)] = random_simple(19)

        mask = vcat(falses(2), trues(5), falses(12))
        @test_throws DimensionMismatch seq[mask] = random_simple(4)
        @test_throws DimensionMismatch seq[mask] = random_simple(6)

        test_bool_arr(seq, mask)
        for i in 1:5
            seq = random_simple(20)
            mask = rand(Bool, 20)
            test_bool_arr(seq, mask)
        end
    end

    @testset "Setindex with colon" begin
        seq1 = SimpleSeq(map(RNA, collect("AGAUGCUCUUAGAC")))
        seq2 = SimpleSeq(map(RNA, collect("AGUCGUAUAUAGGC")))

        seq3 = copy(seq1)
        seq3[:] = seq2
        @test seq3 == seq2

        seq3 = copy(seq2)
        seq3[:] = seq1
        @test seq3 == seq1

        seq3 = copy(seq2)
        seq3[:] = collect("AGAUGCUCUUAGAC")
        @test seq3 == seq1

        seq3 = empty(seq1)[:]
        seq3[:] = RNA[]
        @test isempty(seq3)
    end

    @testset "Setindex with unit range" begin
        seq = SimpleSeq(map(RNA, collect("AUGCUGUAUUCGGAAA"))) # 16 bp

        @test_throws BoundsError seq[15:17] = [RNA_A, RNA_C, RNA_U]
        @test_throws BoundsError seq[5:25] = ""

        @test_throws DimensionMismatch seq[5:7] = "AUCG"

        seq2 = copy(seq)
        seq2[5:4] = ""
        @test seq2 == seq
        
        seq[4:6] = "CAU"
        @test seq == SimpleSeq(map(RNA, collect("AUGCAUUAUUCGGAAA")))
    end

    @testset "Setindex with integer array" begin
        seq = SimpleSeq(map(RNA, collect("AUGCUGCGUAUGUUCUU"))) # 17 bp

        @test_throws BoundsError seq[[0, 1, 2]] = "UUU"
        @test_throws BoundsError seq[[6, 1, 18]] = "UUU"

        @test_throws DimensionMismatch seq[[4,5,6]] = "UU"
        @test_throws DimensionMismatch seq[[4,5,6]] = "ACGU"

        seq[[2,6,8]] = "ACU"
        @test seq == SimpleSeq(map(RNA, collect("AAGCUCCUUAUGUUCUU")))

        seq[1:2:17] = "ACG"^3
        @test seq == SimpleSeq(map(RNA, collect("AACCGCAUCAGGAUCUG")))
    end
end

@testset "BitIndex" begin
    ind = BioSequences.BitIndex{4, UInt64}(16)
    BioSequences.BitsPerSymbol(ind) = BioSequences.BitsPerSymbol{4}()
    BioSequences.bitwidth(UInt64) = 64
    BioSequences.bitwidth(UInt16) = 16
    BioSequences.prevposition(ind) == BioSequences.BitIndex{4, UInt64}(12)
    BioSequences.nextposition(ind) == BioSequences.BitIndex{4, UInt64}(20)
end