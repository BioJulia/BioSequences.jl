@testset "Translation" begin
    # crummy string translation to test against
    standard_genetic_code_dict = Dict{String,Char}(
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

        # translatable ambiguities in the standard code
        "CUN" => 'L', "CCN" => 'P', "CGN" => 'R', "ACN" => 'T',
        "GUN" => 'V', "GCN" => 'A', "GGN" => 'G', "UCN" => 'S'
    )

    function string_translate(seq::AbstractString)
        @assert length(seq) % 3 == 0
        aaseq = Vector{Char}(undef, div(length(seq), 3))
        for i in 1:3:length(seq) - 3 + 1
            aaseq[div(i, 3) + 1] = standard_genetic_code_dict[seq[i:i+2]]
        end
        return String(aaseq)
    end

    function check_translate(seq::AbstractString)
        return string_translate(seq) == String(translate(RNASequence(seq)))
    end

    global reps = 10
    for len in [1, 10, 32, 1000, 10000, 100000]
        @test all(Bool[check_translate(random_translatable_rna(len)) for _ in 1:reps])
    end

    # ambiguous codons
    @test translate(rna"YUGMGG") == aa"LR"
    @test translate(rna"GAYGARGAM") == aa"DEX"

    # BioSequences{RNAAlphabet{2}}
    @test translate(GeneralSequence{RNAAlphabet{2}}("AAAUUUGGGCCC")) == translate(rna"AAAUUUGGGCCC")
    @test translate(GeneralSequence{DNAAlphabet{2}}("AAATTTGGGCCC")) == translate(dna"AAATTTGGGCCC")

    # DNASequence
    @test translate(dna"ATGTAA") == aa"M*"

    # Alternative start codons
    @test translate(rna"GUGUAA", alternative_start = true) == aa"M*"

    @test_throws Exception translate(rna"ACGUACGU")  # can't translate non-multiples of three
    # can't translate N
    @test_throws Exception translate(rna"ACGUACGNU", allow_ambiguous_codons=false)

    # issue #133
    @test translate(rna"GAN") == aa"X"
end
