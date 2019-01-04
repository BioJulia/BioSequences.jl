global reps = 10

@testset "Construction and Conversions" begin
    @test DNACodon(DNA_A, DNA_G, DNA_T) === DNAKmer("AGT")
    @test RNACodon(RNA_A, RNA_G, RNA_U) === RNAKmer("AGU")

    # Check that skipmers and kmers in uint64s survive round trip conversion:
    #   UInt64 → Skipmer / Kmer → UInt64
    function check_uint64_convertion(T::Type, n::UInt64, len::Int)
        return UInt64(Skipmer{T, 2, 3, len}(n)) === n && UInt64(Kmer{T, len}(n)) === n
    end
    
    # Check that skipmers and kmers in uint128's survive round trip conversion:
    #   UInt128 → BigSkipmer / BigKmer → UInt128
    function check_uint128_convertion(T::Type, n::UInt128, len::Int)
        return UInt128(BigSkipmer{T, 2, 3, len}(n)) === n && UInt128(BigKmer{T, len}(n)) === n
    end

    # Check that kmers in strings survive round trip conversion:
    #   String → Kmer → String
    function check_string_construction(T::Type, seq::AbstractString)
        return String(T(seq)) == uppercase(seq)
    end

    # Check that dnakmers can be constructed from a DNASequence
    #   DNASequence → Kmer → DNASequence
    function check_dnasequence_construction(T::Type, seq::DNASequence)
        return DNASequence(T(seq)) == seq
    end

    # Check that rnakmers can be constructed from a RNASequence
    #   RNASequence → Kmer → RNASequence
    function check_rnasequence_construction(T::Type, seq::RNASequence)
        return RNASequence(T(seq)) == seq
    end

    # Check that kmers can be constructed from a BioSequence
    #   BioSequence → Kmer → BioSequence
    function check_biosequence_construction(T::Type, seq::GeneralSequence)
        return GeneralSequence(T(seq)) == seq
    end

    # Check that kmers can be constructed from an array of nucleotides
    #   Vector{NucleicAcid} → Kmer → Vector{NucleicAcid}
    function check_nucarray_kmer(t::Type, seq::Vector{T}) where T <: NucleicAcid
        return String([convert(Char, c) for c in seq]) == String(t(seq...))
    end

    # Check that kmers in strings survive round trip conversion:
    #   String → BioSequence → Kmer → BioSequence → String
    function check_roundabout_construction(T::Type, A, seq::AbstractString)
        #T = eltype(A)
        return String(GeneralSequence{A}(T(GeneralSequence{A}(seq)))) == uppercase(seq)
    end

    @testset "Skipmer and Kmer conversion" begin
        for len in [1, 16, 32]
            # UInt64 conversions
            @test all(Bool[check_uint64_convertion(DNAAlphabet{2}, rand(UInt64(0):UInt64(UInt64(1) << 2len - 1)), len) for _ in 1:reps])
            @test all(Bool[check_uint64_convertion(RNAAlphabet{2}, rand(UInt64(0):UInt64(UInt64(1) << 2len - 1)), len) for _ in 1:reps])

            # String construction
            @test all(Bool[check_string_construction(Kmer{DNAAlphabet{2}}, random_dna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_string_construction(Kmer{RNAAlphabet{2}}, random_rna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_string_construction(Skipmer{DNAAlphabet{2}}, random_dna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_string_construction(Skipmer{RNAAlphabet{2}}, random_rna_kmer(len)) for _ in 1:reps])

            # DNA/RNASequence Constructions
            @test all(Bool[check_dnasequence_construction(Kmer, DNASequence(random_dna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_rnasequence_construction(Kmer, RNASequence(random_rna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_dnasequence_construction(Skipmer, DNASequence(random_dna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_rnasequence_construction(Skipmer, RNASequence(random_rna_kmer(len))) for _ in 1:reps])

            # BioSequence Construction
            @test all(Bool[check_biosequence_construction(Kmer, DNASequence(random_dna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_biosequence_construction(Kmer, RNASequence(random_rna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_biosequence_construction(Skipmer, DNASequence(random_dna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_biosequence_construction(Skipmer, RNASequence(random_rna_kmer(len))) for _ in 1:reps])

            # Construction from nucleotide arrays
            @test all(Bool[check_nucarray_kmer(Kmer, random_dna_kmer_nucleotides(len)) for _ in 1:reps])
            @test all(Bool[check_nucarray_kmer(Kmer, random_rna_kmer_nucleotides(len)) for _ in 1:reps])
            @test all(Bool[check_nucarray_kmer(Skipmer, random_dna_kmer_nucleotides(len)) for _ in 1:reps])
            @test all(Bool[check_nucarray_kmer(Skipmer, random_rna_kmer_nucleotides(len)) for _ in 1:reps])

            # Roundabout conversions
            @test all(Bool[check_roundabout_construction(Kmer, DNAAlphabet{2}, random_dna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_roundabout_construction(Kmer, DNAAlphabet{4}, random_dna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_roundabout_construction(Kmer, RNAAlphabet{2}, random_rna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_roundabout_construction(Kmer, RNAAlphabet{4}, random_rna_kmer(len)) for _ in 1:reps])
            
            @test all(Bool[check_roundabout_construction(Skipmer, DNAAlphabet{2}, random_dna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_roundabout_construction(Skipmer, DNAAlphabet{4}, random_dna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_roundabout_construction(Skipmer, RNAAlphabet{2}, random_rna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_roundabout_construction(Skipmer, RNAAlphabet{4}, random_rna_kmer(len)) for _ in 1:reps])
        end
    end
    
    @testset "BigSkipmer and BigKmer conversion" begin
        for len in [33, 41, 56, 64]
            # UInt64 conversions
            @test all(Bool[check_uint128_convertion(DNAAlphabet{2}, rand(UInt128(0):UInt128(UInt128(1) << 2len - 1)), len) for _ in 1:reps])
            @test all(Bool[check_uint128_convertion(RNAAlphabet{2}, rand(UInt128(0):UInt128(UInt128(1) << 2len - 1)), len) for _ in 1:reps])

            # String construction
            @test all(Bool[check_string_construction(BigKmer{DNAAlphabet{2}}, random_dna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_string_construction(BigKmer{RNAAlphabet{2}}, random_rna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_string_construction(BigSkipmer{DNAAlphabet{2}}, random_dna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_string_construction(BigSkipmer{RNAAlphabet{2}}, random_rna_kmer(len)) for _ in 1:reps])

            # DNA/RNASequence Constructions
            @test all(Bool[check_dnasequence_construction(BigKmer, DNASequence(random_dna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_rnasequence_construction(BigKmer, RNASequence(random_rna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_dnasequence_construction(BigSkipmer, DNASequence(random_dna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_rnasequence_construction(BigSkipmer, RNASequence(random_rna_kmer(len))) for _ in 1:reps])

            # BioSequence Construction
            @test all(Bool[check_biosequence_construction(BigKmer, DNASequence(random_dna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_biosequence_construction(BigKmer, RNASequence(random_rna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_biosequence_construction(BigSkipmer, DNASequence(random_dna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_biosequence_construction(BigSkipmer, RNASequence(random_rna_kmer(len))) for _ in 1:reps])

            # Construction from nucleotide arrays
            @test all(Bool[check_nucarray_kmer(BigKmer, random_dna_kmer_nucleotides(len)) for _ in 1:reps])
            @test all(Bool[check_nucarray_kmer(BigKmer, random_rna_kmer_nucleotides(len)) for _ in 1:reps])
            @test all(Bool[check_nucarray_kmer(BigSkipmer, random_dna_kmer_nucleotides(len)) for _ in 1:reps])
            @test all(Bool[check_nucarray_kmer(BigSkipmer, random_rna_kmer_nucleotides(len)) for _ in 1:reps])

            # Roundabout conversions
            @test all(Bool[check_roundabout_construction(BigKmer, DNAAlphabet{2}, random_dna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_roundabout_construction(BigKmer, DNAAlphabet{4}, random_dna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_roundabout_construction(BigKmer, RNAAlphabet{2}, random_rna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_roundabout_construction(BigKmer, RNAAlphabet{4}, random_rna_kmer(len)) for _ in 1:reps])
            
            @test all(Bool[check_roundabout_construction(BigSkipmer, DNAAlphabet{2}, random_dna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_roundabout_construction(BigSkipmer, DNAAlphabet{4}, random_dna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_roundabout_construction(BigSkipmer, RNAAlphabet{2}, random_rna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_roundabout_construction(BigSkipmer, RNAAlphabet{4}, random_rna_kmer(len)) for _ in 1:reps])
        end
    end

    @test_throws MethodError Kmer() # can't construct 0-mer using `Kmer()`
    @test_throws ArgumentError Kmer(dna"") # 0-mers not allowed
    @test_throws ArgumentError DNAKmer{0}(UInt64(0)) # 0-mers not allowed
    @test_throws ArgumentError RNAKmer{0}(UInt64(0)) # 0-mers not allowed
    @test_throws ArgumentError Kmer(RNA_A, RNA_C, RNA_G, RNA_N, RNA_U) # no Ns in kmers
    @test_throws ArgumentError Kmer(DNA_A, DNA_C, DNA_G, DNA_N, DNA_T) # no Ns in kmers
    @test_throws ArgumentError Kmer(rna"ACGNU")# no Ns in kmers
    @test_throws ArgumentError RNAKmer(rna"ACGNU")# no Ns in kmers
    @test_throws ArgumentError Kmer(dna"ACGNT") # no Ns in kmers
    @test_throws ArgumentError DNAKmer(dna"ACGNT") # no Ns in kmers
    @test_throws MethodError Kmer(RNA_A, DNA_A) # no mixing of RNA and DNA
    @test_throws ArgumentError Kmer{RNAAlphabet{2}}(random_rna(33)) # no kmer larger than 32nt
    @test_throws ArgumentError Kmer{DNAAlphabet{2}}(random_dna(33)) # no kmer larger than 32nt
    @test_throws ArgumentError Kmer(
                      RNA_A, RNA_C, RNA_G, RNA_U, # no kmer larger than 32nt
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U)
    @test_throws ArgumentError Kmer(
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
