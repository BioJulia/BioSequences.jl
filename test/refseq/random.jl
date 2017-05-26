@testset "Random sequence" begin
    for _ in 1:10
        @test randdnaseq(0) == dna""
        @test randrnaseq(0) == rna""
        @test randaaseq(0) == aa""
    end

    for _ in 1:100
        len = rand(0:1000)
        @test length(randdnaseq(len)) == len
        @test length(randrnaseq(len)) == len
        @test length(randaaseq(len)) == len
    end

    for _ in 1:100
        len = 10_000

        counts = countmap(collect(randdnaseq(len)))
        counts_sum = 0
        for nt in dna"ACGT"
            # cdf(Binomial(10_000, 0.25), 2300) ≈ 1.68e-6
            @test 2300 < counts[nt]
            counts_sum += counts[nt]
        end
        @test counts_sum == len

        counts = countmap(collect(randrnaseq(len)))
        counts_sum = 0
        for nt in rna"ACGU"
            # cdf(Binomial(10_000, 0.25), 2300) ≈ 1.68e-6
            @test 2300 < counts[nt]
            counts_sum += counts[nt]
        end
        @test counts_sum == len

        counts = countmap(collect(randaaseq(len)))
        counts_sum = 0
        for aa in aa"ARNDCQEGHILKMFPSTWYV"
            # cdf(Binomial(10000, 0.05), 400) ≈ 1.21e-6
            @test 400 < counts[aa]
            counts_sum += counts[aa]
        end
        @test counts_sum == len
    end
end
