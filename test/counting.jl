@testset "Counting" begin

    @testset "GC content" begin
        @test gc_content(dna"") === 0.0
        @test gc_content(dna"AATA") === 0.0
        @test gc_content(dna"ACGT") === 0.5
        @test gc_content(dna"CGGC") === 1.0
        @test gc_content(dna"ACATTGTGTATAACAAAAGG") === 6 / 20
        @test gc_content(dna"GAGGCGTTTATCATC"[2:end]) === 6 / 14
        
        @test gc_content(rna"") === 0.0
        @test gc_content(rna"AAUA") === 0.0
        @test gc_content(rna"ACGU") === 0.5
        @test gc_content(rna"CGGC") === 1.0
        @test gc_content(rna"ACAUUGUGUAUAACAAAAGG") === 6 / 20
        @test gc_content(rna"GAGGCGUUUAUCAUC"[2:end]) === 6 / 14
        
        @test_throws Exception gc_content(aa"ARN")
        
        Random.seed!(1234)
        for _ in 1:200
            s = randdnaseq(rand(1:100))
            @test gc_content(s) === count(isGC, s) / length(s)
            @test gc_content(LongSequence{DNAAlphabet{2}}(s)) === count(isGC, s) / length(s)
            i = rand(1:lastindex(s))
            j = rand(i-1:lastindex(s))
            @test gc_content(s[i:j]) === (j < i ? 0.0 : count(isGC, s[i:j]) / (j - i + 1))
        end
    end
    
    function testcounter(
        pred::Function,
        alias::Function,
        seqa::BioSequence,
        seqb::BioSequence,
        singlearg::Bool
    )
        # Test that order does not matter.
        @test count(pred, seqa, seqb) == count(pred, seqb, seqa)
        @test BioSequences.count_naive(pred, seqa, seqb) == BioSequences.count_naive(pred, seqb, seqa)
        @test alias(seqa, seqb) == alias(seqb, seqa)
        # Test that result is the same as counting naively.
        @test count(pred, seqa, seqb) == BioSequences.count_naive(pred, seqa, seqb)
        @test count(pred, seqb, seqa) == BioSequences.count_naive(pred, seqb, seqa)
        # Test that the alias function works.
        @test count(pred, seqa, seqb) == alias(seqa, seqb)
        @test count(pred, seqb, seqa) == alias(seqb, seqa)
        if singlearg
            @test count(pred, seqa) == count(pred, (i for i in seqa))
            @test BioSequences.count_naive(pred, seqa) == BioSequences.count(pred, seqa)
        end
    end
    
    function counter_random_tests(
        pred::Function,
        alias::Function,
        alphx::Type{<:Alphabet},
        alphy::Type{<:Alphabet},
        subset::Bool,
        singlearg::Bool
    )
        for _ in 1:10
            seqA = random_seq(alphx, rand(10:100))
            seqB = random_seq(alphy, rand(10:100))
            sa = seqA
            sb = seqB
            if subset
                intA = random_interval(1, length(seqA))
                intB = random_interval(1, length(seqB))
                subA = view(seqA, intA)
                subB = view(seqB, intB)
                sa = subA
                sb = subB
            end
            testcounter(pred, alias, sa, sb, singlearg)
        end
    end
    
    @testset "Mismatches" begin
        for a in (DNAAlphabet, RNAAlphabet)
            # Can't promote views
            for sub in (true, false)
                for n in (4, 2)
                    counter_random_tests(!=, mismatches, a{n}, a{n}, sub, false)
                end
            end
            counter_random_tests(!=, mismatches, a{4}, a{2}, false, false)
            counter_random_tests(!=, mismatches, a{2}, a{4}, false, false)
        end
    end
    
    @testset "Matches" begin
        for a in (DNAAlphabet, RNAAlphabet)
            for sub in (true, false)
                for n in (4, 2)
                    counter_random_tests(==, matches, a{n}, a{n}, sub, false)
                end
            end
            counter_random_tests(==, matches, a{4}, a{2}, false, false)
            counter_random_tests(==, matches, a{2}, a{4}, false, false)
        end
    end
    
    @testset "Ambiguous" begin
        for a in (DNAAlphabet, RNAAlphabet)
            # Can't promote views
            for n in (4, 2)
                for sub in (true, false)
                    counter_random_tests(isambiguous, n_ambiguous, a{n}, a{n}, sub, true)
                end
            end
            counter_random_tests(isambiguous, n_ambiguous, a{4}, a{2}, false, true)
            counter_random_tests(isambiguous, n_ambiguous, a{2}, a{4}, false, true)
        end
    end
    
    @testset "Certain" begin
        for a in (DNAAlphabet, RNAAlphabet)
            for n in (4, 2)
                for sub in (true, false)
                    counter_random_tests(iscertain, n_certain, a{n}, a{n}, sub, true)
                end
            end
            counter_random_tests(iscertain, n_certain, a{4}, a{2}, false, true)
            counter_random_tests(iscertain, n_certain, a{2}, a{4}, false, true)
        end
    end
    
    @testset "Gap" begin
        for a in (DNAAlphabet, RNAAlphabet)
            for n in (4, 2)
                for sub in (true, false)
                    counter_random_tests(isgap, n_gaps, a{n}, a{n}, sub, true)
                end
            end
            counter_random_tests(isgap, n_gaps, a{4}, a{2}, false, true)
            counter_random_tests(isgap, n_gaps, a{2}, a{4}, false, true)
        end
    end
    
end
