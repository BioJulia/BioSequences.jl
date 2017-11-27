@testset "Ambiguous nucleotide positions" begin
    function test_ns(A, seq)
        expected = Int[]
        for i in 1:length(seq)
            if seq[i] == 'N'
                push!(expected, i)
            end
        end
        bioseq = MutableBioSequence{A}(seq)
        @test collect(ambiguous_positions(bioseq)) == expected
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
    for (i, p) in enumerate(ambiguous_positions(dna_seq))
        @test p == pos[i]
    end

    dna_seq = dna"NATTCGRATY"
    pos = [1, 7, 10]
    for (i, p) in enumerate(ambiguous_positions(dna_seq))
        @test p == pos[i]
    end

    rna_seq = rna"ANANANA"
    pos = [2, 4, 6]
    for (i, p) in enumerate(ambiguous_positions(rna_seq))
        @test p == pos[i]
    end

    rna_seq = rna"NAUUCGRAUY"
    pos = [1, 7, 10]
    for (i, p) in enumerate(ambiguous_positions(rna_seq))
        @test p == pos[i]
    end
end
