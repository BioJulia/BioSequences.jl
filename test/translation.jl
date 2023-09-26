@testset "Translation" begin
    # crummy string translation to test against
    standard_genetic_code_dict = let
        d_ = Dict{String,Char}(
            "AAA" => 'K', "AAC" => 'N', "AAG" => 'K', "AAU" => 'N',
            "ACA" => 'T', "ACC" => 'T', "ACG" => 'T', "ACU" => 'T',
            "AGA" => 'R', "AGC" => 'S', "AGG" => 'R', "AGU" => 'S',
            "AUA" => 'I', "AUC" => 'I', "AUG" => 'M', "AUU" => 'I',
            "CAA" => 'Q', "CAC" => 'H', "CAG" => 'Q', "CAU" => 'H',
            "CCA" => 'P', "CCC" => 'P', "CCG" => 'P', "CCU" => 'P',
            "CGA" => 'R', "CGC" => 'R', "CGG" => 'R', "CGU" => 'R',
            "CUA" => 'L', "CUC" => 'L', "CUG" => 'L', "CUU" => 'L',
            "GAA" => 'E', "GAC" => 'D', "GAG" => 'E', "GAU" => 'D',
            "GCA" => 'A', "GCC" => 'A', "GCG" => 'A', "GCU" => 'A',
            "GGA" => 'G', "GGC" => 'G', "GGG" => 'G', "GGU" => 'G',
            "GUA" => 'V', "GUC" => 'V', "GUG" => 'V', "GUU" => 'V',
            "UAA" => '*', "UAC" => 'Y', "UAG" => '*', "UAU" => 'Y',
            "UCA" => 'S', "UCC" => 'S', "UCG" => 'S', "UCU" => 'S',
            "UGA" => '*', "UGC" => 'C', "UGG" => 'W', "UGU" => 'C',
            "UUA" => 'L', "UUC" => 'F', "UUG" => 'L', "UUU" => 'F',
        )
        d = Dict{NTuple{3, RNA}, Char}(map(RNA, Tuple(k)) => v for (k, v) in d_)

        for (a, b, c) in Iterators.product(ntuple(i -> alphabet(RNA), Val{3}())...)
            (a == RNA_Gap || b == RNA_Gap || c == RNA_Gap) && continue
            (a, b, c) ∈ keys(d) && continue
            possible = ' '
            for (x, y, z) in Iterators.product(map(BioSequences.UnambiguousRNAs, (a, b, c))...)
                new_possible = d[(x, y, z)]
                if possible == ' '
                    possible = new_possible
                elseif possible == new_possible
                    nothing
                elseif possible ∈ ('J', 'I', 'L') && new_possible ∈ ('I', 'L')
                    possible = 'J'
                elseif possible ∈ ('B', 'D', 'N') && new_possible ∈ ('D', 'N')
                    possible = 'B'
                elseif possible ∈ ('Z', 'Q', 'E') && new_possible ∈ ('Q', 'E')
                    possible = 'Z'
                else
                    possible = 'X'
                    break
                end
            end
            d[(a, b, c)] = possible
        end
        result = Dict(join(k) => v for (k, v) in d)
        for (k, v) in result
            if 'U' in k
                result[replace(k, 'U'=>'T')] = v
            end
        end
        result
    end

    function string_translate(seq::AbstractString)
        @assert length(seq) % 3 == 0
        aaseq = Vector{Char}(undef, div(length(seq), 3))
        for i in 1:3:length(seq) - 3 + 1
            aaseq[div(i, 3) + 1] = standard_genetic_code_dict[seq[i:i+2]]
        end
        return String(aaseq)
    end

    # Basics
    @test length(BioSequences.standard_genetic_code) == 64
    buf = IOBuffer()
    show(buf, BioSequences.standard_genetic_code)
    @test !iszero(length(take!(buf))) # just test it doesn't error

    # TransTables
    @test BioSequences.TransTables() isa BioSequences.TransTables
    @test BioSequences.ncbi_trans_table[1] === BioSequences.standard_genetic_code
    buf = IOBuffer()
    show(buf, BioSequences.ncbi_trans_table)
    @test !iszero(length(take!(buf)))

    # Test translation values
    sampler = SamplerWeighted(rna"ACGUMRSVWYHKDBN", [fill(0.2, 4); fill(0.018181818, 10)])
    for len in [3, 9, 54, 102]
        for A in [DNAAlphabet{4}(), RNAAlphabet{4}()]
            for attempt in 1:25
                seq = randseq(A, sampler, len)
                @test string(translate(seq)) == string_translate(string(seq))
                for window in [1:3, 5:13, 11:67]
                    checkbounds(Bool, eachindex(seq), window) || continue
                    v = view(seq, window)
                    @test string(translate(v)) == string_translate(string(v))
                end
            end
        end
        for A in [DNAAlphabet{2}(), RNAAlphabet{2}()]
            seq = randseq(A, len)
            @test string(translate(seq)) == string_translate(string(seq))
        end
    end

    # ambiguous codons
    @test translate(rna"YUGMGG") == aa"LR"
    @test translate(rna"GAYGARGAM") == aa"DEX"
    @test translate(rna"MUCGGG") == aa"JG"
    @test translate(rna"AAASAAUUU") == aa"KZF"

    # BioSequences{RNAAlphabet{2}}
    @test translate(LongSequence{RNAAlphabet{2}}("AAAUUUGGGCCC")) == translate(rna"AAAUUUGGGCCC")
    @test translate(LongSequence{DNAAlphabet{2}}("AAATTTGGGCCC")) == translate(dna"AAATTTGGGCCC")

    # LongDNA{4}
    @test translate(dna"ATGTAA") == aa"M*"

    # Alternative start codons
    @test translate(rna"GUGUAA", alternative_start = true) == aa"M*"

    @test_throws Exception translate(rna"ACGUACGU")  # can't translate non-multiples of three
    # can't translate N
    @test_throws Exception translate(rna"ACGUACGNU", allow_ambiguous_codons=false)

    # Can't translate gaps
    @test_throws Exception translate(dna"A-G")
    @test_throws Exception translate(dna"---")
    @test_throws Exception translate(dna"AACGAT-A-")

    # issue #133
    @test translate(rna"GAN") == aa"X"

    # Test views
    seq = dna"TAGTGCNTAGDACGGGWAAABCCGTAC"
    @test translate(view(seq, 1:0)) == aa""
    @test translate(LongSubSeq{RNAAlphabet{4}}(seq, 2:length(seq)-2)) == translate(seq[2:end-2])
end
