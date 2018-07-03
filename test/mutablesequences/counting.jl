@testset "Site counting" begin

    dna_alphabets = (DNAAlphabet{4}, DNAAlphabet{2})
    rna_alphabets = (RNAAlphabet{4}, RNAAlphabet{2})

    @testset "Specific count methods" begin

        function generate_possibilities_tester{A<:NucleicAcidAlphabet}(::Type{A}) where A<:NucleicAcidAlphabet
            # TODO: Resolve this use of symbols(A()).
            symbolset = symbols(A())
            arra = Vector{eltype(A)}()
            arrb = Vector{eltype(A)}()
            for i in 1:length(symbolset), j in i:length(symbolset)
                push!(arra, symbolset[i])
                push!(arrb, symbolset[j])
            end
            return GeneralSequence{A}(arra), GeneralSequence{A}(arrb)
        end

        for alphset in (dna_alphabets, rna_alphabets)
            # Answers to these tests were worked out manually to verify
            # counting is working correctly.
            # seqA and seqB contain all possible observations of sites.

            # 4 bit encoded sequences
            seqA, seqB = generate_possibilities_tester(alphset[1])
            @test count(Certain, seqA, seqB) == count(Certain, seqB, seqA) == 10
            @test count(Gap, seqA, seqB) == count(Gap, seqB, seqA) == 16
            @test count(Ambiguous, seqA, seqB) == count(Ambiguous, seqB, seqA) == 121
            @test count(Match, seqA, seqB) == count(Match, seqB, seqA) == length(symbols(alphset[1]()))
            @test count(Mismatch, seqA, seqB) == count(Mismatch, seqB, seqA) == (length(seqA) - length(symbols(alphset[1]())))

            # 2 bit encoded sequences
            seqA, seqB = generate_possibilities_tester(alphset[2])
            @test count(Certain, seqA, seqB) == count(Certain, seqB, seqA) == 10
            @test count(Gap, seqA, seqB) == count(Gap, seqB, seqA) == 0
            @test count(Ambiguous, seqA, seqB) == count(Ambiguous, seqB, seqA) == 0
            @test count(Match, seqA, seqB) == count(Match, seqB, seqA) == length(symbols(alphset[2]()))
            @test count(Mismatch, seqA, seqB) == count(Mismatch, seqB, seqA) == (length(seqA) - length(symbols(alphset[2]())))
        end
    end

    @testset "Randomized tests" begin

        # A test counting function which is naive.
        @inline function testcount(::Type{P}, a::GeneralSequence, b::GeneralSequence) where {P<:BioSequences.Position}
            k = 0
            @inbounds for idx in 1:min(lastindex(a), lastindex(b))
                k += issite(P, a, b, idx)
            end
            return k
        end
        issite(::Type{Ambiguous}, a::GeneralSequence, idx) = isambiguous(a[idx])
        @inline function issite(::Type{Ambiguous}, a::GeneralSequence, b::GeneralSequence, idx)
            return issite(Ambiguous, a, idx) | issite(Ambiguous, b, idx)
        end
        issite(::Type{Certain}, a::GeneralSequence, idx) = iscertain(a[idx])
        @inline function issite(::Type{Certain}, a::GeneralSequence, b::GeneralSequence, idx)
            return issite(Certain, a, idx) & issite(Certain, b, idx)
        end
        issite(::Type{Gap}, a::GeneralSequence, idx) = isgap(a[idx])
        @inline function issite(::Type{Gap}, a::GeneralSequence, b::GeneralSequence, idx)
            return issite(Gap, a, idx) | issite(Gap, b, idx)
        end
        issite(::Type{Match}, a::GeneralSequence, b::GeneralSequence, idx) = a[idx] == b[idx]
        issite(::Type{Mismatch}, a::GeneralSequence, b::GeneralSequence, idx) = a[idx] != b[idx]

        # Randomized tests get performed with a naive counting function
        # which is intuitive and works, but that is nowhere near as quick.

        function testcounting(::Type{S}, a, b) where S<:Site
            @test count(S, a, b) == count(S, b, a) == testcount(S, a, b)
        end

        function testforencs(a::Int, b::Int, subset::Bool)
            for alphabet in (DNAAlphabet, RNAAlphabet)
                for _ in  1:50
                    #println("TESTING SEQUENCES:")
                    seqA = random_seq(alphabet{a}, rand(10:100))
                    seqB = random_seq(alphabet{b}, rand(10:100))
                    #=
                    println("A seq: ", seqA)
                    println("B seq: ", seqB)
                    println("A length: ", length(seqA))
                    println("B length: ", length(seqB))
                    =#
                    sa = seqA
                    sb = seqB
                    if subset
                        intA = random_interval(1, length(seqA))
                        intB = random_interval(1, length(seqB))
                        subA = seqA[intA]
                        subB = seqB[intB]
                        #=
                        println("A subset: ", intA)
                        println("B subset: ", intB)
                        println("A subseq: ", subA)
                        println("B subseq: ", subB)
                        =#
                        sa = subA
                        sb = subB
                    end
                    #println("Testing Mismatch")
                    testcounting(Mismatch, sa, sb)
                    #println("Testing Match")
                    testcounting(Match, sa, sb)
                    #println("Testing Certain")
                    testcounting(Certain, sa, sb)
                    #println("Testing Gap")
                    testcounting(Gap, sa, sb)
                    #println("Testing Ambiguous")
                    testcounting(Ambiguous, sa, sb)
                end
            end
        end

        @testset "Testing 4-bit seq pairs" begin
            @testset "Full random sequences" begin
                testforencs(4, 4, false)
            end
            @testset "Subset random sequences" begin
                testforencs(4, 4, true)
            end
        end

        @testset "Testing 2-bit seq pairs" begin
            @testset "Full random sequences" begin
                testforencs(2, 2, false)
            end
            @testset "Subset random sequences" begin
                testforencs(2, 2, true)
            end
        end

        @testset "Testing mixed bit seq pairs" begin
            @testset "Full random sequences" begin
                testforencs(4, 2, false)
                testforencs(2, 4, false)
            end
            @testset "Subset random sequences" begin
                testforencs(4, 2, true)
                testforencs(2, 4, true)
            end
        end
    end

    @testset "Pairwise methods" begin
        @testset "4-bit encoded sequences" begin
            dnas = [dna"ATCGCCA-", dna"ATCGCCTA", dna"ATCGCCT-", dna"GTCGCCTA"]
            rnas = [rna"AUCGCCA-", rna"AUCGCCUA", rna"AUCGCCU-", rna"GUCGCCUA"]
            answer_mismatch = [0 2 1 3;
                               2 0 1 1;
                               1 1 0 2;
                               3 1 2 0]
            answer_match = [0 6 7 5;
                            6 0 7 7;
                            7 7 0 6;
                            5 7 6 0]
            for i in (dnas, rnas)
                @test count_pairwise(Mismatch, i...) == answer_mismatch
                @test count_pairwise(Match, i...) == answer_match
                @test count_pairwise(Certain, i...) == [0 7 7 7;
                                                        7 0 7 8;
                                                        7 7 0 7;
                                                        7 8 7 0]
                @test count_pairwise(Ambiguous, i...) == [0 0 0 0;
                                                          0 0 0 0;
                                                          0 0 0 0;
                                                          0 0 0 0]
                @test count_pairwise(Gap, i...) == [0 1 1 1;
                                                    1 0 1 0;
                                                    1 1 0 1;
                                                    1 0 1 0]
            end
        end
        @testset "2-bit encoded sequences" begin
            dnas = [GeneralSequence{DNAAlphabet{2}}("ATCGCCAC"),
                    GeneralSequence{DNAAlphabet{2}}("ATCGCCTA"),
                    GeneralSequence{DNAAlphabet{2}}("ATCGCCTT"),
                    GeneralSequence{DNAAlphabet{2}}("GTCGCCTA")]
            rnas = [GeneralSequence{RNAAlphabet{2}}("AUCGCCAC"),
                    GeneralSequence{RNAAlphabet{2}}("AUCGCCUA"),
                    GeneralSequence{RNAAlphabet{2}}("AUCGCCUU"),
                    GeneralSequence{RNAAlphabet{2}}("GUCGCCUA")]
            answer_mismatch = [0 2 2 3;
                               2 0 1 1;
                               2 1 0 2;
                               3 1 2 0]
            answer_match = [0 6 6 5;
                            6 0 7 7;
                            6 7 0 6;
                            5 7 6 0]
            for i in (dnas, rnas)
                @test count_pairwise(Mismatch, i...) == answer_mismatch
                @test count_pairwise(Match, i...) == answer_match
                @test count_pairwise(Certain, i...) == [0 8 8 8;
                                                        8 0 8 8;
                                                        8 8 0 8;
                                                        8 8 8 0]
                @test count_pairwise(Ambiguous, i...) == [0 0 0 0;
                                                          0 0 0 0;
                                                          0 0 0 0;
                                                          0 0 0 0]
                @test count_pairwise(Gap, i...) == [0 0 0 0;
                                                    0 0 0 0;
                                                    0 0 0 0;
                                                    0 0 0 0]
            end
        end
    end

    @testset "Windowed methods" begin
        @testset "4-bit encoded sequences" begin
            dnaA = dna"ATCGCCA-M"
            dnaB = dna"ATCGCCTAA"
            rnaA = rna"AUCGCCA-M"
            rnaB = rna"AUCGCCUAA"

            for seqs in ((dnaA, dnaB), (rnaA, rnaB))
                @test count(Certain, seqs[1], seqs[2], 3, 1) == [IntervalValue(1, 3, 3),
                                                                 IntervalValue(2, 4, 3),
                                                                 IntervalValue(3, 5, 3),
                                                                 IntervalValue(4, 6, 3),
                                                                 IntervalValue(5, 7, 3),
                                                                 IntervalValue(6, 8, 2),
                                                                 IntervalValue(7, 9, 1)]
                @test count(Ambiguous, seqs[1], seqs[2], 3, 1) == [IntervalValue(1, 3, 0),
                                                                   IntervalValue(2, 4, 0),
                                                                   IntervalValue(3, 5, 0),
                                                                   IntervalValue(4, 6, 0),
                                                                   IntervalValue(5, 7, 0),
                                                                   IntervalValue(6, 8, 0),
                                                                   IntervalValue(7, 9, 1)]
                @test count(Gap, seqs[1], seqs[2], 3, 1) == [IntervalValue(1, 3, 0),
                                                             IntervalValue(2, 4, 0),
                                                             IntervalValue(3, 5, 0),
                                                             IntervalValue(4, 6, 0),
                                                             IntervalValue(5, 7, 0),
                                                             IntervalValue(6, 8, 1),
                                                             IntervalValue(7, 9, 1)]
                @test count(Match, seqs[1], seqs[2], 3, 1) == [IntervalValue(1, 3, 3),
                                                               IntervalValue(2, 4, 3),
                                                               IntervalValue(3, 5, 3),
                                                               IntervalValue(4, 6, 3),
                                                               IntervalValue(5, 7, 2),
                                                               IntervalValue(6, 8, 1),
                                                               IntervalValue(7, 9, 0)]
                @test count(Mismatch, seqs[1], seqs[2], 3, 1) == [IntervalValue(1, 3, 0),
                                                                  IntervalValue(2, 4, 0),
                                                                  IntervalValue(3, 5, 0),
                                                                  IntervalValue(4, 6, 0),
                                                                  IntervalValue(5, 7, 1),
                                                                  IntervalValue(6, 8, 2),
                                                                  IntervalValue(7, 9, 3)]
            end
        end
        @testset "2-bit encoded sequences" begin
            dnaA = GeneralSequence{DNAAlphabet{2}}("ATCGCCATT")
            dnaB = GeneralSequence{DNAAlphabet{2}}("ATCGCCTAA")
            rnaA = GeneralSequence{RNAAlphabet{2}}("AUCGCCAUU")
            rnaB = GeneralSequence{RNAAlphabet{2}}("AUCGCCUAA")

            for seqs in ((dnaA, dnaB), (rnaA, rnaB))
                @test count(Certain, seqs[1], seqs[2], 3, 1) == [IntervalValue(1, 3, 3),
                                                                 IntervalValue(2, 4, 3),
                                                                 IntervalValue(3, 5, 3),
                                                                 IntervalValue(4, 6, 3),
                                                                 IntervalValue(5, 7, 3),
                                                                 IntervalValue(6, 8, 3),
                                                                 IntervalValue(7, 9, 3)]
                @test count(Ambiguous, seqs[1], seqs[2], 3, 1) == [IntervalValue(1, 3, 0),
                                                                   IntervalValue(2, 4, 0),
                                                                   IntervalValue(3, 5, 0),
                                                                   IntervalValue(4, 6, 0),
                                                                   IntervalValue(5, 7, 0),
                                                                   IntervalValue(6, 8, 0),
                                                                   IntervalValue(7, 9, 0)]
                @test count(Gap, seqs[1], seqs[2], 3, 1) == [IntervalValue(1, 3, 0),
                                                             IntervalValue(2, 4, 0),
                                                             IntervalValue(3, 5, 0),
                                                             IntervalValue(4, 6, 0),
                                                             IntervalValue(5, 7, 0),
                                                             IntervalValue(6, 8, 0),
                                                             IntervalValue(7, 9, 0)]
                @test count(Match, seqs[1], seqs[2], 3, 1) == [IntervalValue(1, 3, 3),
                                                               IntervalValue(2, 4, 3),
                                                               IntervalValue(3, 5, 3),
                                                               IntervalValue(4, 6, 3),
                                                               IntervalValue(5, 7, 2),
                                                               IntervalValue(6, 8, 1),
                                                               IntervalValue(7, 9, 0)]
                @test count(Mismatch, seqs[1], seqs[2], 3, 1) == [IntervalValue(1, 3, 0),
                                                                  IntervalValue(2, 4, 0),
                                                                  IntervalValue(3, 5, 0),
                                                                  IntervalValue(4, 6, 0),
                                                                  IntervalValue(5, 7, 1),
                                                                  IntervalValue(6, 8, 2),
                                                                  IntervalValue(7, 9, 3)]
            end
        end
        @testset "Mixed encodings" begin
            dnaA = dna"ATCGCCA-M"
            dnaB = GeneralSequence{DNAAlphabet{2}}("ATCGCCTAA")
            rnaA = rna"AUCGCCA-M"
            rnaB = GeneralSequence{RNAAlphabet{2}}("AUCGCCUAA")

            for seqs in ((dnaA, dnaB), (rnaA, rnaB))
                @test count(Certain, seqs[1], seqs[2], 3, 1) == [IntervalValue(1, 3, 3),
                                                                 IntervalValue(2, 4, 3),
                                                                 IntervalValue(3, 5, 3),
                                                                 IntervalValue(4, 6, 3),
                                                                 IntervalValue(5, 7, 3),
                                                                 IntervalValue(6, 8, 2),
                                                                 IntervalValue(7, 9, 1)]
                @test count(Ambiguous, seqs[1], seqs[2], 3, 1) == [IntervalValue(1, 3, 0),
                                                                   IntervalValue(2, 4, 0),
                                                                   IntervalValue(3, 5, 0),
                                                                   IntervalValue(4, 6, 0),
                                                                   IntervalValue(5, 7, 0),
                                                                   IntervalValue(6, 8, 0),
                                                                   IntervalValue(7, 9, 1)]
                @test count(Gap, seqs[1], seqs[2], 3, 1) == [IntervalValue(1, 3, 0),
                                                             IntervalValue(2, 4, 0),
                                                             IntervalValue(3, 5, 0),
                                                             IntervalValue(4, 6, 0),
                                                             IntervalValue(5, 7, 0),
                                                             IntervalValue(6, 8, 1),
                                                             IntervalValue(7, 9, 1)]
                @test count(Match, seqs[1], seqs[2], 3, 1) == [IntervalValue(1, 3, 3),
                                                               IntervalValue(2, 4, 3),
                                                               IntervalValue(3, 5, 3),
                                                               IntervalValue(4, 6, 3),
                                                               IntervalValue(5, 7, 2),
                                                               IntervalValue(6, 8, 1),
                                                               IntervalValue(7, 9, 0)]
                @test count(Mismatch, seqs[1], seqs[2], 3, 1) == [IntervalValue(1, 3, 0),
                                                                  IntervalValue(2, 4, 0),
                                                                  IntervalValue(3, 5, 0),
                                                                  IntervalValue(4, 6, 0),
                                                                  IntervalValue(5, 7, 1),
                                                                  IntervalValue(6, 8, 2),
                                                                  IntervalValue(7, 9, 3)]
            end
        end
    end
end
