global reps = 10

@testset "Construction and Conversions" begin
    @test DNACodon(DNA_A, DNA_G, DNA_T) === DNAKmer("AGT")
    @test RNACodon(RNA_A, RNA_G, RNA_U) === RNAKmer("AGU")

    # Check that kmers in strings survive round trip conversion:
    #   UInt64 → Kmer → UInt64
    function check_uint64_convertion(T::Type, n::UInt64, len::Int)
        return UInt64(Kmer{T, len}(n)) === n
    end

    # Check that kmers in strings survive round trip conversion:
    #   String → Kmer → String
    function check_string_construction(T::Type, seq::AbstractString)
        return String(Kmer{T}(seq)) == uppercase(seq)
    end

    # Check that dnakmers can be constructed from a DNASequence
    #   DNASequence → Kmer → DNASequence
    function check_dnasequence_construction(seq::DNASequence)
        return DNASequence(DNAKmer(seq)) == seq
    end

    # Check that rnakmers can be constructed from a RNASequence
    #   RNASequence → Kmer → RNASequence
    function check_rnasequence_construction(seq::RNASequence)
        return RNASequence(RNAKmer(seq)) == seq
    end

    # Check that kmers can be constructed from a BioSequence
    #   BioSequence → Kmer → BioSequence
    function check_biosequence_construction(seq::GeneralSequence)
        return GeneralSequence(Kmer(seq)) == seq
    end

    # Check that kmers can be constructed from an array of nucleotides
    #   Vector{NucleicAcid} → Kmer → Vector{NucleicAcid}
    function check_nucarray_kmer(seq::Vector{T}) where T <: NucleicAcid
        return String([convert(Char, c) for c in seq]) == String(Kmer(seq...))
    end

    # Check that kmers in strings survive round trip conversion:
    #   String → BioSequence → Kmer → BioSequence → String
    function check_roundabout_construction(A, seq::AbstractString)
        T = eltype(A)
        return String(GeneralSequence{A}(Kmer(GeneralSequence{A}(seq)))) == uppercase(seq)
    end

    for len in [1, 16, 32]
        # UInt64 conversions
        @test all(Bool[check_uint64_convertion(DNA, rand(UInt64(0):UInt64(UInt64(1) << 2len - 1)), len) for _ in 1:reps])
        @test all(Bool[check_uint64_convertion(RNA, rand(UInt64(0):UInt64(UInt64(1) << 2len - 1)), len) for _ in 1:reps])

        # String construction
        @test all(Bool[check_string_construction(DNA, random_dna_kmer(len)) for _ in 1:reps])
        @test all(Bool[check_string_construction(RNA, random_rna_kmer(len)) for _ in 1:reps])

        # DNA/RNASequence Constructions
        @test all(Bool[check_dnasequence_construction(DNASequence(random_dna_kmer(len))) for _ in 1:reps])
        @test all(Bool[check_rnasequence_construction(RNASequence(random_rna_kmer(len))) for _ in 1:reps])

        # BioSequence Construction
        @test all(Bool[check_biosequence_construction(DNASequence(random_dna_kmer(len))) for _ in 1:reps])
        @test all(Bool[check_biosequence_construction(RNASequence(random_rna_kmer(len))) for _ in 1:reps])

        # Construction from nucleotide arrays
        @test all(Bool[check_nucarray_kmer(random_dna_kmer_nucleotides(len)) for _ in 1:reps])
        @test all(Bool[check_nucarray_kmer(random_rna_kmer_nucleotides(len)) for _ in 1:reps])

        # Roundabout conversions
        @test all(Bool[check_roundabout_construction(DNAAlphabet{2}, random_dna_kmer(len)) for _ in 1:reps])
        @test all(Bool[check_roundabout_construction(DNAAlphabet{4}, random_dna_kmer(len)) for _ in 1:reps])
        @test all(Bool[check_roundabout_construction(RNAAlphabet{2}, random_rna_kmer(len)) for _ in 1:reps])
        @test all(Bool[check_roundabout_construction(RNAAlphabet{4}, random_rna_kmer(len)) for _ in 1:reps])
    end

    @test_throws Exception Kmer() # can't construct 0-mer using `Kmer()`
    @test_throws Exception Kmer(dna"") # 0-mers not allowed
    @test_throws Exception DNAKmer{0}(UInt64(0)) # 0-mers not allowed
    @test_throws Exception RNAKmer{0}(UInt64(0)) # 0-mers not allowed
    @test_throws Exception Kmer(RNA_A, RNA_C, RNA_G, RNA_N, RNA_U) # no Ns in kmers
    @test_throws Exception Kmer(DNA_A, DNA_C, DNA_G, DNA_N, DNA_T) # no Ns in kmers
    @test_throws Exception Kmer(rna"ACGNU")# no Ns in kmers
    @test_throws Exception RNAKmer(rna"ACGNU")# no Ns in kmers
    @test_throws Exception Kmer(dna"ACGNT") # no Ns in kmers
    @test_throws Exception DNAKmer(dna"ACGNT") # no Ns in kmers
    @test_throws Exception Kmer(RNA_A, DNA_A) # no mixing of RNA and DNA
    @test_throws Exception Kmer(random_rna(33)) # no kmer larger than 32nt
    @test_throws Exception Kmer(random_dna(33)) # no kmer larger than 32nt
    @test_throws Exception Kmer(
                      RNA_A, RNA_C, RNA_G, RNA_U, # no kmer larger than 32nt
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U)
    @test_throws Exception Kmer(
                      DNA_A, DNA_C, DNA_G, DNA_T, # no kmer larger than 32nt
                      DNA_A, DNA_C, DNA_G, DNA_T,
                      DNA_A, DNA_C, DNA_G, DNA_T,
                      DNA_A, DNA_C, DNA_G, DNA_T,
                      DNA_A, DNA_C, DNA_G, DNA_T,
                      DNA_A, DNA_C, DNA_G, DNA_T,
                      DNA_A, DNA_C, DNA_G, DNA_T,
                      DNA_A, DNA_C, DNA_G, DNA_T,
                      DNA_A, DNA_C, DNA_G, DNA_T)

    @testset "From strings" begin
        @test DNAKmer("ACTG") == Kmer(DNASequence("ACTG"))
        @test RNAKmer("ACUG") == Kmer(RNASequence("ACUG"))

        # N is not allowed in Kmers
        @test_throws Exception DNAKmer("ACGTNACGT")
        @test_throws Exception RNAKmer("ACGUNACGU")

        # Test string literals
        @test kmer"ACTG" == Kmer(DNASequence("ACTG"))
        @test isa(kmer"ACGT", DNAKmer{4})
        @test_throws LoadError eval(:(kmer"ACGN"))
        @test_throws LoadError eval(:(kmer"ACG-"))
    end
    
    @testset "Capacity" begin
        @test BioSequences.capacity(DNAKmer(random_dna_kmer(10))) == 32
        @test BioSequences.capacity(RNAKmer(random_rna_kmer(10))) == 32
        @test BioSequences.capacity(DNAKmer(random_dna_kmer(32))) == 32
        @test BioSequences.capacity(RNAKmer(random_rna_kmer(32))) == 32
    end
    
    @testset "N unused" begin
        @test BioSequences.n_unused(DNAKmer(random_dna_kmer(10))) == 22
        @test BioSequences.n_unused(RNAKmer(random_rna_kmer(10))) == 22
        @test BioSequences.n_unused(DNAKmer(random_dna_kmer(32))) == 0
        @test BioSequences.n_unused(RNAKmer(random_rna_kmer(32))) == 0
    end
    
    @testset "Span" begin
        @test BioSequences.span(DNAKmer(random_dna_kmer(10))) == 10
        @test BioSequences.span(RNAKmer(random_rna_kmer(10))) == 10
        @test BioSequences.span(DNAKmer(random_dna_kmer(32))) == 32
        @test BioSequences.span(RNAKmer(random_rna_kmer(32))) == 32
    end
end
