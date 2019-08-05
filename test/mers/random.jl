@testset "Random" begin
    @testset for k in 1:64
        if k <= 32
            for _ in 1:10
                kmer = rand(DNAMer{k})
                @test isa(kmer, DNAMer{k})
                kmer = rand(RNAMer{k})
                @test isa(kmer, RNAMer{k})
            end
            for size in [0, 1, 2, 5, 10, 100]
                @test length(rand(DNAMer{k}, size)) == size
                @test length(rand(RNAMer{k}, size)) == size
            end
            kmers = rand(DNAMer{k}, 10_000)
            for i in 1:k
                a = sum([kmer[i] for kmer in kmers] .== DNA_A)
                c = sum([kmer[i] for kmer in kmers] .== DNA_C)
                g = sum([kmer[i] for kmer in kmers] .== DNA_G)
                t = sum([kmer[i] for kmer in kmers] .== DNA_T)
                @test 2200 ≤ a ≤ 2800
                @test 2200 ≤ c ≤ 2800
                @test 2200 ≤ g ≤ 2800
                @test 2200 ≤ t ≤ 2800
                @test a + c + g + t == 10_000
            end
        end
        
        for _ in 1:10
            kmer = rand(BigDNAMer{k})
            @test isa(kmer, BigDNAMer{k})
            kmer = rand(BigRNAMer{k})
            @test isa(kmer, BigRNAMer{k})
        end
        for size in [0, 1, 2, 5, 10, 100]
            @test length(rand(BigDNAMer{k}, size)) == size
            @test length(rand(BigRNAMer{k}, size)) == size
        end
        kmers = rand(BigDNAMer{k}, 10_000)
        for i in 1:k
            a = sum([kmer[i] for kmer in kmers] .== DNA_A)
            c = sum([kmer[i] for kmer in kmers] .== DNA_C)
            g = sum([kmer[i] for kmer in kmers] .== DNA_G)
            t = sum([kmer[i] for kmer in kmers] .== DNA_T)
            @test 2200 ≤ a ≤ 2800
            @test 2200 ≤ c ≤ 2800
            @test 2200 ≤ g ≤ 2800
            @test 2200 ≤ t ≤ 2800
            @test a + c + g + t == 10_000
        end
    end
end
