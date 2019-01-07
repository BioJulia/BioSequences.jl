global reps = 10

@testset "Construction and Conversions" begin
    @test DNACodon(DNA_A, DNA_G, DNA_T) === DNAKmer("AGT")
    @test RNACodon(RNA_A, RNA_G, RNA_U) === RNAKmer("AGU")
    
    # Check that kmers in strings survive round trip conversion:
    #   String → Kmer → String
    function check_string_construction(::Type{Skipmer{U, A, M, N, K}}, seq::AbstractString) where {U, A, M, N, K}
        return String(Skipmer{U, A, M, N, K}(seq)) == uppercase(seq)
    end

    # Check that dnakmers can be constructed from a DNASequence
    #   DNASequence → Kmer → DNASequence
    function check_dnasequence_construction(::Type{Skipmer{U, A, M, N, K}}, seq::DNASequence) where {U, A, M, N, K}
        return DNASequence(Skipmer{U, A, M, N, K}(seq)) == seq
    end

    # Check that rnakmers can be constructed from a RNASequence
    #   RNASequence → Kmer → RNASequence
    function check_rnasequence_construction(::Type{Skipmer{U, A, M, N, K}}, seq::RNASequence) where {U, A, M, N, K}
        return RNASequence(Skipmer{U, A, M, N, K}(seq)) == seq
    end

    # Check that kmers can be constructed from a BioSequence
    #   BioSequence → Kmer → BioSequence
    function check_biosequence_construction(::Type{Skipmer{U, A, M, N, K}}, seq::LongSequence) where {U, A, M, N, K}
        return LongSequence(Skipmer{U, A, M, N, K}(seq)) == seq
    end

    # Check that kmers can be constructed from an array of nucleotides
    #   Vector{NucleicAcid} → Kmer → Vector{NucleicAcid}
    function check_nucarray_kmer(::Type{Skipmer{U, A, M, N, K}}, seq::Vector{T}) where {T <: NucleicAcid, U, A, M, N, K}
        return String([convert(Char, c) for c in seq]) == String(Skipmer{U, A, M, N, K}(seq))
    end

    # Check that kmers in strings survive round trip conversion:
    #   String → BioSequence → Kmer → BioSequence → String
    function check_roundabout_construction(::Type{Skipmer{U, A, M, N, K}}, A2, seq::AbstractString) where {U, A, M, N, K}
        return String(LongSequence{A2}(Skipmer{U, A, M, N, K}(LongSequence{A2}(seq)))) == uppercase(seq)
    end
    
    function check_uint_conversion(::Type{Skipmer{U, A, M, N, K}}) where {U, A, M, N, K}
        uint = rand(typemin(U):U(one(U) << 2K - 1))
        return convert(U, Skipmer{U, A, M, N, K}(uint)) === uint
    end

    @testset "Skipmer conversion" begin
        for len in [1, 16, 32]
            # UInt64 conversions
            @test all(Bool[check_uint_conversion(Skipmer{UInt64, DNAAlphabet{2}, 2, 3, len}) for _ in 1:reps])
            @test all(Bool[check_uint_conversion(Skipmer{UInt128, DNAAlphabet{2}, 2, 3, len}) for _ in 1:reps])
            @test all(Bool[check_uint_conversion(Skipmer{UInt64, RNAAlphabet{2}, 2, 3, len}) for _ in 1:reps])
            @test all(Bool[check_uint_conversion(Skipmer{UInt128, RNAAlphabet{2}, 2, 3, len}) for _ in 1:reps])

            # String construction
            @test all(Bool[check_string_construction(DNAKmer{len}, random_dna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_string_construction(RNAKmer{len}, random_rna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_string_construction(Skipmer{UInt64, DNAAlphabet{2}, 2, 3, len}, random_dna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_string_construction(Skipmer{UInt64, RNAAlphabet{2}, 2, 3, len}, random_rna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_string_construction(Skipmer{UInt128, DNAAlphabet{2}, 2, 3, len}, random_dna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_string_construction(Skipmer{UInt128, RNAAlphabet{2}, 2, 3, len}, random_rna_kmer(len)) for _ in 1:reps])

            # DNA/RNASequence Constructions
            @test all(Bool[check_dnasequence_construction(DNAKmer{len}, DNASequence(random_dna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_rnasequence_construction(RNAKmer{len}, RNASequence(random_rna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_dnasequence_construction(Skipmer{UInt64, DNAAlphabet{2}, 2, 3, len}, DNASequence(random_dna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_rnasequence_construction(Skipmer{UInt64, RNAAlphabet{2}, 2, 3, len}, RNASequence(random_rna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_dnasequence_construction(Skipmer{UInt128, DNAAlphabet{2}, 2, 3, len}, DNASequence(random_dna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_rnasequence_construction(Skipmer{UInt128, RNAAlphabet{2}, 2, 3, len}, RNASequence(random_rna_kmer(len))) for _ in 1:reps])

            # BioSequence Construction
            @test all(Bool[check_biosequence_construction(DNAKmer{len}, DNASequence(random_dna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_biosequence_construction(RNAKmer{len}, RNASequence(random_rna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_biosequence_construction(Skipmer{UInt64, DNAAlphabet{2}, 2, 3, len}, DNASequence(random_dna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_biosequence_construction(Skipmer{UInt64, RNAAlphabet{2}, 2, 3, len}, RNASequence(random_rna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_biosequence_construction(Skipmer{UInt128, DNAAlphabet{2}, 2, 3, len}, DNASequence(random_dna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_biosequence_construction(Skipmer{UInt128, RNAAlphabet{2}, 2, 3, len}, RNASequence(random_rna_kmer(len))) for _ in 1:reps])

            # Construction from nucleotide arrays
            @test all(Bool[check_nucarray_kmer(DNAKmer{len}, random_dna_kmer_nucleotides(len)) for _ in 1:reps])
            @test all(Bool[check_nucarray_kmer(RNAKmer{len}, random_rna_kmer_nucleotides(len)) for _ in 1:reps])
            @test all(Bool[check_nucarray_kmer(Skipmer{UInt64, DNAAlphabet{2}, 2, 3, len}, random_dna_kmer_nucleotides(len)) for _ in 1:reps])
            @test all(Bool[check_nucarray_kmer(Skipmer{UInt64, RNAAlphabet{2}, 2, 3, len}, random_rna_kmer_nucleotides(len)) for _ in 1:reps])
            @test all(Bool[check_nucarray_kmer(Skipmer{UInt128, DNAAlphabet{2}, 2, 3, len}, random_dna_kmer_nucleotides(len)) for _ in 1:reps])
            @test all(Bool[check_nucarray_kmer(Skipmer{UInt128, RNAAlphabet{2}, 2, 3, len}, random_rna_kmer_nucleotides(len)) for _ in 1:reps])

            # Roundabout conversions
            @test all(Bool[check_roundabout_construction(DNAKmer{len}, DNAAlphabet{2}, random_dna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_roundabout_construction(DNAKmer{len}, DNAAlphabet{4}, random_dna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_roundabout_construction(RNAKmer{len}, RNAAlphabet{2}, random_rna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_roundabout_construction(RNAKmer{len}, RNAAlphabet{4}, random_rna_kmer(len)) for _ in 1:reps])
            
            @test all(Bool[check_roundabout_construction(Skipmer{UInt64, DNAAlphabet{2}, 2, 3, len}, DNAAlphabet{2}, random_dna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_roundabout_construction(Skipmer{UInt64, DNAAlphabet{2}, 2, 3, len}, DNAAlphabet{4}, random_dna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_roundabout_construction(Skipmer{UInt64, RNAAlphabet{2}, 2, 3, len}, RNAAlphabet{2}, random_rna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_roundabout_construction(Skipmer{UInt64, RNAAlphabet{2}, 2, 3, len}, RNAAlphabet{4}, random_rna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_roundabout_construction(Skipmer{UInt128, DNAAlphabet{2}, 2, 3, len}, DNAAlphabet{2}, random_dna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_roundabout_construction(Skipmer{UInt128, DNAAlphabet{2}, 2, 3, len}, DNAAlphabet{4}, random_dna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_roundabout_construction(Skipmer{UInt128, RNAAlphabet{2}, 2, 3, len}, RNAAlphabet{2}, random_rna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_roundabout_construction(Skipmer{UInt128, RNAAlphabet{2}, 2, 3, len}, RNAAlphabet{4}, random_rna_kmer(len)) for _ in 1:reps])
        end
    end

    @test_throws MethodError Kmer() # can't construct 0-mer using `Kmer()`
    @test_throws ArgumentError DNAKmer(dna"") # 0-mers not allowed
    @test_throws ArgumentError DNAKmer{0}(UInt64(0)) # 0-mers not allowed
    @test_throws ArgumentError RNAKmer{0}(UInt64(0)) # 0-mers not allowed
    @test_throws ArgumentError RNAKmer(RNA_A, RNA_C, RNA_G, RNA_N, RNA_U) # no Ns in kmers
    @test_throws ArgumentError DNAKmer(DNA_A, DNA_C, DNA_G, DNA_N, DNA_T) # no Ns in kmers
    @test_throws ArgumentError RNAKmer(rna"ACGNU")# no Ns in kmers
    @test_throws ArgumentError DNAKmer(dna"ACGNT") # no Ns in kmers
    @test_throws MethodError DNAKmer(RNA_A, DNA_A) # no mixing of RNA and DNA
    @test_throws ArgumentError Kmer{UInt64, RNAAlphabet{2}}(random_rna(33)) # no kmer larger than 32nt
    @test_throws ArgumentError Kmer{UInt64, DNAAlphabet{2}}(random_dna(33)) # no kmer larger than 32nt
    @test_throws ArgumentError Kmer{UInt64, RNAAlphabet{2}}(
                      RNA_A, RNA_C, RNA_G, RNA_U, # no kmer larger than 32nt
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U)
    @test_throws ArgumentError Kmer{UInt64, DNAAlphabet{2}}(
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
        @test DNAKmer("ACTG") == DNAKmer(DNASequence("ACTG"))
        @test RNAKmer("ACUG") == RNAKmer(RNASequence("ACUG"))

        # N is not allowed in Kmers
        @test_throws Exception DNAKmer("ACGTNACGT")
        @test_throws Exception RNAKmer("ACGUNACGU")

        # Test string literals
        @test kmer"ACTG" == DNAKmer(DNASequence("ACTG"))
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
