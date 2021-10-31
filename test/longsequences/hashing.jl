@testset "Hash" begin
    s = dna"ACGTACGT"
    v = view(s, 1:lastindex(s))

	for seq in [s, v]
	    @test isa(hash(seq), UInt64)
	    @test hash(seq) === hash(dna"ACGTACGT")
	    @test hash(seq) !== hash(seq[1:6])
	    @test hash(seq) !== hash(seq[1:7])
	    @test hash(seq) === hash(seq[1:8])
	    @test hash(seq[1:4]) === hash(dna"ACGT")
	    @test hash(seq[2:4]) === hash(dna"CGT")
	    @test hash(seq[3:4]) === hash(dna"GT")
	    @test hash(seq[4:4]) === hash(dna"T")
	    @test hash(seq[5:8]) === hash(dna"ACGT")
	end

	@test hash(s) == hash(v)
	@test hash(s[2:4]) == hash(v[2:4])

	for i in 1:4
		s1 = LongDNA{4}("A"^(i-1))
		s2 = LongDNA{4}("A"^i)
		@test hash(s1) != hash(s2)
		v1 = view(s2, 1:lastindex(s2) - 1)
		v2 = view(s2, 1:lastindex(s2))
		@test hash(v1) != hash(v2)
	end

    for n in [1, 2, 3, 6, 8, 9, 11, 15, 16, 20], seq in [dna"A", dna"AC", dna"ACG", dna"ACGT", dna"ACGTN"]
        @test hash(seq^n) === hash((dna""     * seq^n)[1:end])
        @test hash(seq^n) === hash((dna"T"    * seq^n)[2:end])
        @test hash(seq^n) === hash((dna"TT"   * seq^n)[3:end])
        @test hash(seq^n) === hash((dna"TTT"  * seq^n)[4:end])
        @test hash(seq^n) === hash((dna"TTTT" * seq^n)[5:end])

        @test hash(seq^n) === hash((seq^n * dna""    )[1:end  ])
        @test hash(seq^n) === hash((seq^n * dna"T"   )[1:end-1])
        @test hash(seq^n) === hash((seq^n * dna"TT"  )[1:end-2])
        @test hash(seq^n) === hash((seq^n * dna"TTT" )[1:end-3])
        @test hash(seq^n) === hash((seq^n * dna"TTTT")[1:end-4])
    end

    @test hash(rna"AAUU") === hash(rna"AAUU")
    @test hash(rna"AAUUAA"[3:5]) === hash(rna"UUA")
    @test hash(aa"MTTQAPMFTQPLQ") === hash(aa"MTTQAPMFTQPLQ")
    @test hash(aa"MTTQAPMFTQPLQ"[5:10]) === hash(aa"APMFTQ")

	# Test hash of longer view to engange some inner loops
	seq = randdnaseq(250)
	@test hash(seq[33:201]) == hash(view(seq, 33:201))
end
