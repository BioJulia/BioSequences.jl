@testset "Condition based iteration" begin
    @testset "Ambiguous nucleotide positions" begin
        function test_ns(A, seq)
            expected = Tuple{Int, eltype(A)}[]
            for i in 1:length(seq)
                if seq[i] == 'N'
                    push!(expected, (i, DNA_N))
                end
            end
            bioseq = LongSequence{A}(seq)
            @test collect(each(isambiguous, bioseq)) == expected
        end

        for len in [1, 10, 32, 1000, 10000, 100000]
            test_ns(DNAAlphabet{4}, random_dna(len))
            test_ns(RNAAlphabet{4}, random_rna(len))

            probs = [0.25, 0.25, 0.25, 0.25, 0.00]
            test_ns(DNAAlphabet{2}, random_dna(len, probs))
            test_ns(RNAAlphabet{2}, random_rna(len, probs))
        end

        dna_seq = dna"ANANANA"
        pos = [2, 4, 6]
        j = 1
        for (p, nuc) in each(isambiguous, dna_seq)
            @test p == pos[j]
            j += 1
        end

        dna_seq = dna"NATTCGRATY"
        pos = [1, 7, 10]
        j = 1
        for (p, nuc) in each(isambiguous, dna_seq)
            @test p == pos[j]
            j += 1
        end

        rna_seq = rna"ANANANA"
        pos = [2, 4, 6]
        j = 1
        for (p, nuc) in each(isambiguous, rna_seq)
            @test p == pos[j]
            j += 1
        end

        rna_seq = rna"NAUUCGRAUY"
        pos = [1, 7, 10]
        j = 1
        for (p, nuc) in each(isambiguous,rna_seq)
            @test p == pos[j]
            j += 1
        end

        @test collect(each(isambiguous, rna_seq)) == collect(zip(pos, (RNA(i) for i in "NRY")))
    end
end
