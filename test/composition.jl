@testset "Composition" begin
    function string_nucleotide_count(::Type{DNA}, seq::AbstractString)
        counts = Dict{DNA, Int}(
            DNA_A => 0,
            DNA_C => 0,
            DNA_G => 0,
            DNA_T => 0,
            DNA_N => 0 )
        for c in seq
            counts[convert(DNA, c)] += 1
        end
        return counts
    end

    function string_nucleotide_count(::Type{RNA}, seq::AbstractString)
        counts = Dict{RNA, Int}(
            RNA_A => 0,
            RNA_C => 0,
            RNA_G => 0,
            RNA_U => 0,
            RNA_N => 0 )
        for c in seq
            counts[convert(RNA, c)] += 1
        end
        return counts
    end

    function check_nucleotide_count(::Type{DNA}, seq::AbstractString)
        string_counts = string_nucleotide_count(DNA, seq)
        seq_counts = composition(DNASequence(seq))
        return string_counts[DNA_A] == seq_counts[DNA_A] &&
               string_counts[DNA_C] == seq_counts[DNA_C] &&
               string_counts[DNA_G] == seq_counts[DNA_G] &&
               string_counts[DNA_T] == seq_counts[DNA_T] &&
               string_counts[DNA_N] == seq_counts[DNA_N]
    end

    function check_nucleotide_count(::Type{RNA}, seq::AbstractString)
        string_counts = string_nucleotide_count(RNA, seq)
        seq_counts = composition(RNASequence(seq))
        return string_counts[RNA_A] == seq_counts[RNA_A] &&
               string_counts[RNA_C] == seq_counts[RNA_C] &&
               string_counts[RNA_G] == seq_counts[RNA_G] &&
               string_counts[RNA_U] == seq_counts[RNA_U] &&
               string_counts[RNA_N] == seq_counts[RNA_N]
    end

    function check_kmer_nucleotide_count(::Type{DNA}, seq::AbstractString)
        string_counts = string_nucleotide_count(DNA, seq)
        kmer_counts = composition(DNAKmer(seq))
        return string_counts[DNA_A] == kmer_counts[DNA_A] &&
               string_counts[DNA_C] == kmer_counts[DNA_C] &&
               string_counts[DNA_G] == kmer_counts[DNA_G] &&
               string_counts[DNA_T] == kmer_counts[DNA_T] &&
               string_counts[DNA_N] == kmer_counts[DNA_N]
    end

    function check_kmer_nucleotide_count(::Type{RNA}, seq::AbstractString)
        string_counts = string_nucleotide_count(RNA, seq)
        kmer_counts = composition(RNAKmer(seq))
        return string_counts[RNA_A] == kmer_counts[RNA_A] &&
               string_counts[RNA_C] == kmer_counts[RNA_C] &&
               string_counts[RNA_G] == kmer_counts[RNA_G] &&
               string_counts[RNA_U] == kmer_counts[RNA_U] &&
               string_counts[RNA_N] == kmer_counts[RNA_N]
    end

    global reps = 10
    for len in [1, 10, 32, 1000, 10000, 100000]
        @test all(Bool[check_nucleotide_count(DNA, random_dna(len)) for _ in 1:reps])
        @test all(Bool[check_nucleotide_count(RNA, random_rna(len)) for _ in 1:reps])
    end

    @test composition(aa"MTTQAPMFTQPLQSVVV")[AA_E] === 0
    @test composition(aa"MTTQAPMFTQPLQSVVV")[AA_A] === 1
    @test composition(aa"MTTQAPMFTQPLQSVVV")[AA_P] === 2
    @test composition(aa"MTTQAPMFTQPLQSVVV")[AA_V] === 3

    let comp = composition(ReferenceSequence(dna"AGTTATGN"))
        @test comp[DNA_A] === 2
        @test comp[DNA_C] === 0
        @test comp[DNA_G] === 2
        @test comp[DNA_T] === 3
        @test comp[DNA_N] === 1
    end

    for len in [1, 10, 32]
        @test all(Bool[check_kmer_nucleotide_count(DNA, random_dna_kmer(len)) for _ in 1:reps])
        @test all(Bool[check_kmer_nucleotide_count(RNA, random_rna_kmer(len)) for _ in 1:reps])
    end

    comp = composition(dna"ACGT")
    @test merge!(comp, composition(dna"ACGT")) === comp
    @test comp == composition(dna"ACGT"^2)

    for _ in 1:30
        m, n = rand(0:1000), rand(0:1000)
        seq1 = DNASequence(random_dna(m))
        seq2 = DNASequence(random_dna(n))
        @test composition(seq1 * seq2) == merge(composition(seq1), composition(seq2))
    end

    comp = composition(each(DNAKmer{2}, dna"ACGTACGT"))
    @test comp[DNAKmer("AC")] == 2
    @test comp[DNAKmer("CG")] == 2
    @test comp[DNAKmer("GT")] == 2
    @test comp[DNAKmer("TA")] == 1
    @test comp[DNAKmer("AA")] == 0

    comp = composition(each(DNAKmer{3}, dna"ACGTACGT"))
    @test comp[DNAKmer("ACG")] == 2
    @test comp[DNAKmer("CGT")] == 2
    @test comp[DNAKmer("GTA")] == 1
    @test comp[DNAKmer("TAC")] == 1
    @test comp[DNAKmer("AAA")] == 0

    comp = composition(DNASequence[dna"ATCG", dna"GCTA", dna"ATCGG",
                                   dna"ATCG", dna"ATCG", dna"GCTA"])
    @test comp[dna"ATCG"] == 3
    @test comp[dna"GCTA"] == 2
    @test comp[dna"ATCGG"] == 1
end
