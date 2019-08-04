@testset "Shuffle" begin
    for s in ["A", "C", "G", "T"]
        kmer = DNAMer(s)
        bigkmer = BigDNAMer(s)
        @test kmer === shuffle(kmer)
        @test bigkmer === shuffle(bigkmer)
    end

    function count(kmer)
        a = c = g = t = 0
        for x in kmer
            a += x == DNA_A
            c += x == DNA_C
            g += x == DNA_G
            t += x == DNA_T
        end
        return a, c, g, t
    end

    for k in 1:64, _ in 1:10
        if k <= 32
            kmer = rand(DNAMer{k})
            @test count(kmer) == count(shuffle(kmer))
            if k ≥ 30
                @test kmer != shuffle(kmer)
            end
        end
        kmer = rand(BigDNAMer{k})
        @test count(kmer) == count(shuffle(kmer))
        if k ≥ 30
            @test kmer != shuffle(kmer)
        end
    end
end
