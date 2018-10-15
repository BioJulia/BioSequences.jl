@testset "Random" begin
    @testset for k in 1:32
        for _ in 1:10
            kmer = rand(DNAKmer{k})
            @test isa(kmer, DNAKmer{k})
            kmer = rand(RNAKmer{k})
            @test isa(kmer, RNAKmer{k})
        end

        for size in [0, 1, 2, 5, 10, 100]
            @test length(rand(DNAKmer{k}, size)) == size
            @test length(rand(RNAKmer{k}, size)) == size
        end

        kmers = rand(DNAKmer{k}, 10_000)
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
