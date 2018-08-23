@testset "Demultiplexer" begin
    function randdna(n)
        return DNASequence(rand([DNA_A, DNA_C, DNA_G, DNA_T, DNA_N], n))
    end

    function make_errors(seq, p=0.03)
        seq = copy(seq)
        nucs = DNA['A', 'C', 'G', 'T', 'N']
        i = 1
        while i ≤ lastindex(seq)
            if rand() < p
                r = rand()
                if r < 0.1
                    # insertion
                    insert!(seq, i, rand(nucs))
                elseif r < 0.2
                    # deletion
                    deleteat!(seq, i)
                else
                    # substitution
                    seq[i] = rand(nucs)
                end
            end
            i += 1
        end
        return seq
    end

    @testset "Hamming distance" begin
        barcodes = DNASequence["ATGG", "CAGA", "GGAA", "TACG"]
        dplxr = Demultiplexer(barcodes, n_max_errors=1, distance=:hamming)

        for i in 1:lastindex(barcodes)
            @test dplxr[i] == barcodes[i]
        end

        @test demultiplex(dplxr, dna"ATGG") === (1, 0)
        @test demultiplex(dplxr, dna"CAGA") === (2, 0)
        @test demultiplex(dplxr, dna"GGAA") === (3, 0)
        @test demultiplex(dplxr, dna"TACG") === (4, 0)

        # every 1bp substitution is recoverable
        for (i, barcode) in enumerate(barcodes)
            # substitution
            for j in 1:lastindex(barcode)
                barcode′ = copy(barcode)
                for nt in [DNA_A, DNA_C, DNA_G, DNA_T, DNA_N]
                    barcode′[j] = nt
                    @test demultiplex(dplxr, barcode′) == (i, count(Mismatch, barcode, barcode′))
                end
            end
        end

        n_ok = 0
        n_ok_with_fallback = 0
        for _ in 1:10_000
            i = rand(1:lastindex(barcodes))
            seq = make_errors(barcodes[i] * randdna(10))
            n_ok += demultiplex(dplxr, seq)[1] == i
            n_ok_with_fallback += demultiplex(dplxr, seq, true)[1] == i
        end
        # empirically, n_ok / 10_000 is ~0.985
        # @show n_ok
        @test n_ok / 10_000 > 0.98
        @test n_ok < n_ok_with_fallback
    end

    @testset "Levenshtein distance" begin
        barcodes = DNASequence["ATGG", "CAGA", "GGAA", "TACG"]
        dplxr = Demultiplexer(barcodes, n_max_errors=1, distance=:levenshtein)

        for i in 1:lastindex(barcodes)
            @test dplxr[i] == barcodes[i]
        end

        @test demultiplex(dplxr, dna"ATGG") === (1, 0)
        @test demultiplex(dplxr, dna"CAGA") === (2, 0)
        @test demultiplex(dplxr, dna"GGAA") === (3, 0)
        @test demultiplex(dplxr, dna"TACG") === (4, 0)

        # every 1bp substitution/insertion/deletion is recoverable
        for (i, barcode) in enumerate(barcodes)
            # substitution
            for j in 1:lastindex(barcode)
                barcode′ = copy(barcode)
                for nt in [DNA_A, DNA_C, DNA_G, DNA_T, DNA_N]
                    barcode′[j] = nt
                    i′, d = demultiplex(dplxr, barcode′)
                    @test i′ == i && 0 ≤ d ≤ 1
                end
            end
            # insertion
            for j in 1:lastindex(barcode)
                barcode′ = copy(barcode)
                insert!(barcode′, j, rand([DNA_A, DNA_C, DNA_G, DNA_T, DNA_N]))
                i′, d = demultiplex(dplxr, barcode′)
                @test i′ == i && 0 ≤ d ≤ 1
            end
            # deletion
            for j in 1:lastindex(barcode)
                barcode′ = copy(barcode)
                deleteat!(barcode′, j)
                i′, d = demultiplex(dplxr, barcode′)
                @test i′ == i && 0 ≤ d ≤ 1
            end
        end

        n_ok = 0
        n_ok_with_fallback = 0
        for _ in 1:10_000
            i = rand(1:lastindex(barcodes))
            seq = make_errors(barcodes[i] * randdna(10))
            n_ok += demultiplex(dplxr, seq)[1] == i
            n_ok_with_fallback += demultiplex(dplxr, seq, true)[1] == i
        end
        # empirically, n_ok / 10_000 is ~0.995
        # @show n_ok
        @test n_ok / 10_000 > 0.99
        @test n_ok < n_ok_with_fallback
    end
end
