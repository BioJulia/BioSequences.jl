@testset "Position weight matrix" begin
    @testset "PFM" begin
        m = [
            1   2  3
            4   5  6
            7   8  9
            10 11 12
        ]
        pfm = PFM{DNA}(m)
        @test pfm isa PFM{DNA,Int}
        @test convert(Matrix, pfm) == m
        @test all(pfm[i] == m[i] for i in eachindex(pfm))
        @test all(pfm[i,j] == m[i,j] for i in 1:4, j in 1:3)
        @test all(pfm[ACGT[i],j] == m[i,j] for i in 1:4, j in 1:3)
        @test startswith(sprint(show, pfm), string(summary(pfm), ":\n"))
        @test_throws ArgumentError PFM{DNA}(vcat(m, [13 14 15]))

        # indexing
        pfm[5] = 100
        @test pfm[5] == pfm[DNA_A,2] == 100
        pfm[2,3] = 200
        @test pfm[2,3] == pfm[DNA_C,3] == 200
        pfm[DNA_G,1] = 300
        @test pfm[3,1] == pfm[DNA_G,1] == 300
        @test_throws ArgumentError pfm[DNA_N,1]

        # broadcast
        @test log.(pfm) isa PFM{DNA,Float64}
        @test log.(pfm) == PFM{DNA}(log.(m))
        @test pfm .+ 1 isa PFM{DNA,Int}
        @test pfm .+ 1 == PFM{DNA}(m .+ 1)
        @test pfm .+ 1.1 isa PFM{DNA,Float64}
        @test pfm .+ 1.1 == PFM{DNA}(m .+ 1.1)
        @test pfm .+ [0,1,2,3] isa PFM{DNA,Int}
        @test pfm .+ [0,1,2,3] == PFM{DNA}(m .+ [0,1,2,3])

        set = LongDNA{4}.(split(
        """
        ACG
        ATG
        AGG
        ACC
        """))
        pfm = PFM(set)
        @test pfm == [
            4 0 0  # A
            0 2 1  # C
            0 1 3  # G
            0 1 0  # T
        ]
        @test pfm == PFM(Set(set))
        @test_throws ArgumentError PFM(LongDNA{4}[])
        @test_throws ArgumentError PFM(["foo"])
#        @test_throws ArgumentError PFM([LongDNA{4}("AA"), LongDNA{4}("AA")])
        @test_throws ArgumentError PFM([LongDNA{4}("AA"), LongDNA{4}("AAA")])
    end

    @testset "PWM" begin
        m = [
             1  2  3
             4  5  6
             7  8  9
            10 11 12
        ]
        pfm = PFM{DNA}(m)
        pwm = PWM(pfm)
        raw = log2.((pfm ./ sum(pfm, dims=1)) ./ fill(1/4, 4))
        @test pwm isa PWM{DNA,Float64}
        @test pwm == raw
        @test PWM(pfm, prior=[0.1, 0.4, 0.4, 0.1]) == log2.((pfm ./ sum(pfm, dims=1)) ./ [0.1, 0.4, 0.4, 0.1])
        @test PWM(pfm) == PWM{DNA}(raw)
        @test maxscore(pwm) ≈ sum(maximum(raw, dims=1))
        @test maxscore(PWM{DNA}(zeros(4, 0))) === 0.0
        @test PWM(pfm .+ 0.1) isa PWM{DNA,Float64}  # pseudo count
        @test_throws ArgumentError PWM{DNA}(hcat(m, [13, 14, 15]))
        @test_throws ArgumentError PWM(pfm, prior=normalize([0,1,2,3], 1))
        @test_throws ArgumentError PWM(pfm, prior=normalize([0,1,2,3], 1).+1e-3)
        @test all(pwm[i] === pwm[i] for i in eachindex(pwm))
        @test all(pwm[i,j] === raw[i,j] for i in 1:4, j in 1:3)
        @test all(pwm[ACGT[i],j] === pwm[i,j] for i in 1:4, j in 1:3)
        @test startswith(sprint(show, pwm), string(summary(pwm), ":\n"))
    end

    @testset "findfirst and findlast" begin
        seq = dna"ACGATNATCGCGTANTG"
        data = [
            1.0 0.1 0.2
            0.0 0.2 0.3
            0.1 0.2 0.0
            0.9 0.5 0.2
        ]
        pwm = PWM{DNA}(data)
        @test maxscore(pwm) == 1.8
        @test scoreat(seq, pwm, 1) === 1.2
        @test findfirst(PWMSearchQuery(pwm, 1.0), seq) === 1
        @test findfirst(PWMSearchQuery(pwm, 1.4), seq) === 4
        @test findfirst(PWMSearchQuery(pwm, 1.8), seq) === 7
        @test findfirst(PWMSearchQuery(pwm, 2.0), seq) === nothing
        @test_throws ArgumentError findfirst(PWMSearchQuery(pwm, 1.0), LongRNA{4}(seq))
        @test findlast(PWMSearchQuery(pwm, 1.0), seq) == 14
        @test findlast(PWMSearchQuery(pwm, 1.4), seq) === 7
        @test findlast(PWMSearchQuery(pwm, 1.8), seq) === 7
        @test findlast(PWMSearchQuery(pwm, 2.0), seq) === nothing
        @test_throws ArgumentError findlast(PWMSearchQuery(pwm, 1.0), LongRNA{4}(seq))
    end
end
