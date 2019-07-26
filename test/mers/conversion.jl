global reps = 10

@testset "Construction and Conversions" begin
    @test DNACodon(DNA_A, DNA_G, DNA_T) === DNAMer("AGT")
    @test RNACodon(RNA_A, RNA_G, RNA_U) === RNAMer("AGU")
    
    # Check that kmers in strings survive round trip conversion:
    #   String → Kmer → String
    function check_string_construction(::Type{T}, seq::AbstractString) where {T<:AbstractMer}
        return String(T(seq)) == uppercase(seq)
    end

    # Check that DNAMers can be constructed from a LongDNASeq
    #   LongDNASeq → Kmer → LongDNASeq
    function check_dnasequence_construction(::Type{T}, seq::LongDNASeq) where {T<:AbstractMer}
        return LongDNASeq(T(seq)) == seq
    end

    # Check that RNAMers can be constructed from a LongRNASeq
    #   LongRNASeq → Kmer → LongRNASeq
    function check_rnasequence_construction(::Type{T}, seq::LongRNASeq) where {T<:AbstractMer}
        return LongRNASeq(T(seq)) == seq
    end

    # Check that kmers can be constructed from a BioSequence
    #   BioSequence → Kmer → BioSequence
    function check_biosequence_construction(::Type{T}, seq::LongSequence) where {T<:AbstractMer}
        return LongSequence(T(seq)) == seq
    end

    # Check that kmers can be constructed from an array of nucleotides
    #   Vector{NucleicAcid} → Kmer → Vector{NucleicAcid}
    function check_nucarray_kmer(::Type{M}, seq::Vector{T}) where {T<:NucleicAcid,M<:AbstractMer}
        return String([convert(Char, c) for c in seq]) == String(M(seq))
    end

    # Check that kmers in strings survive round trip conversion:
    #   String → BioSequence → Kmer → BioSequence → String
    function check_roundabout_construction(::Type{T}, A2, seq::AbstractString) where {T<:AbstractMer}
        return String(LongSequence{A2}(T(LongSequence{A2}(seq)))) == uppercase(seq)
    end
    
    function check_uint_conversion(::Type{T}) where {T<:AbstractMer}
        U = BioSequences.encoded_data_type(T)
        uint = rand(typemin(U):U(one(U) << 2BioSequences.ksize(T) - 1))
        return convert(U, T(uint)) === uint
    end

    @testset "Skipmer conversion" begin
        for len in [1, 16, 32, 64]
            
            if len <= 32
                # UInt64 conversions
                @test all(Bool[check_uint_conversion(DNAMer{len}) for _ in 1:reps])
                @test all(Bool[check_uint_conversion(RNAMer{len}) for _ in 1:reps])
                # String construction
                @test all(Bool[check_string_construction(DNAMer{len}, random_dna_kmer(len)) for _ in 1:reps])
                @test all(Bool[check_string_construction(RNAMer{len}, random_rna_kmer(len)) for _ in 1:reps])
                # DNA/LongRNASeq Constructions
                @test all(Bool[check_dnasequence_construction(DNAMer{len}, LongDNASeq(random_dna_kmer(len))) for _ in 1:reps])
                @test all(Bool[check_rnasequence_construction(RNAMer{len}, LongRNASeq(random_rna_kmer(len))) for _ in 1:reps])
                # BioSequence Construction
                @test all(Bool[check_biosequence_construction(DNAMer{len}, LongDNASeq(random_dna_kmer(len))) for _ in 1:reps])
                @test all(Bool[check_biosequence_construction(RNAMer{len}, LongRNASeq(random_rna_kmer(len))) for _ in 1:reps])
                # Construction from nucleotide arrays
                @test all(Bool[check_nucarray_kmer(DNAMer{len}, random_dna_kmer_nucleotides(len)) for _ in 1:reps])
                @test all(Bool[check_nucarray_kmer(RNAMer{len}, random_rna_kmer_nucleotides(len)) for _ in 1:reps])
                # Roundabout conversions
                @test all(Bool[check_roundabout_construction(DNAMer{len}, DNAAlphabet{2}, random_dna_kmer(len)) for _ in 1:reps])
                @test all(Bool[check_roundabout_construction(DNAMer{len}, DNAAlphabet{4}, random_dna_kmer(len)) for _ in 1:reps])
                @test all(Bool[check_roundabout_construction(RNAMer{len}, RNAAlphabet{2}, random_rna_kmer(len)) for _ in 1:reps])
                @test all(Bool[check_roundabout_construction(RNAMer{len}, RNAAlphabet{4}, random_rna_kmer(len)) for _ in 1:reps])
            end
            
            # UInt64 conversions
            @test all(Bool[check_uint_conversion(BigDNAMer{len}) for _ in 1:reps])
            @test all(Bool[check_uint_conversion(BigDNAMer{len}) for _ in 1:reps])
            # String construction
            @test all(Bool[check_string_construction(BigDNAMer{len}, random_dna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_string_construction(BigRNAMer{len}, random_rna_kmer(len)) for _ in 1:reps])
            # DNA/LongRNASeq Constructions
            @test all(Bool[check_dnasequence_construction(BigDNAMer{len}, LongDNASeq(random_dna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_rnasequence_construction(BigRNAMer{len}, LongRNASeq(random_rna_kmer(len))) for _ in 1:reps])
            # BioSequence Construction
            @test all(Bool[check_biosequence_construction(BigDNAMer{len}, LongDNASeq(random_dna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_biosequence_construction(BigRNAMer{len}, LongRNASeq(random_rna_kmer(len))) for _ in 1:reps])
            # Construction from nucleotide arrays
            @test all(Bool[check_biosequence_construction(BigDNAMer{len}, LongDNASeq(random_dna_kmer(len))) for _ in 1:reps])
            @test all(Bool[check_biosequence_construction(BigRNAMer{len}, LongRNASeq(random_rna_kmer(len))) for _ in 1:reps])
            # Roundabout conversions
            @test all(Bool[check_roundabout_construction(BigDNAMer{len}, DNAAlphabet{2}, random_dna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_roundabout_construction(BigDNAMer{len}, DNAAlphabet{4}, random_dna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_roundabout_construction(BigRNAMer{len}, RNAAlphabet{2}, random_rna_kmer(len)) for _ in 1:reps])
            @test all(Bool[check_roundabout_construction(BigRNAMer{len}, RNAAlphabet{4}, random_rna_kmer(len)) for _ in 1:reps])
        end
    end

    @test_throws MethodError Mer() # can't construct 0-mer using `Kmer()`
    @test_throws MethodError BigMer() # can't construct 0-mer using `BigKmer()`
    @test_throws ArgumentError DNAMer(dna"") # 0-mers not allowed
    @test_throws ArgumentError BigDNAMer(dna"") # 0-mers not allowed
    @test_throws ArgumentError DNAMer{0}(UInt64(0)) # 0-mers not allowed
    @test_throws ArgumentError RNAMer{0}(UInt64(0)) # 0-mers not allowed
    @test_throws ArgumentError RNAMer(RNA_A, RNA_C, RNA_G, RNA_N, RNA_U) # no Ns in kmers
    @test_throws ArgumentError DNAMer(DNA_A, DNA_C, DNA_G, DNA_N, DNA_T) # no Ns in kmers
    @test_throws ArgumentError RNAMer(rna"ACGNU")# no Ns in kmers
    @test_throws ArgumentError DNAMer(dna"ACGNT") # no Ns in kmers
    @test_throws MethodError DNAMer(RNA_A, DNA_A) # no mixing of RNA and DNA
    @test_throws ArgumentError RNAMer(random_rna(33)) # no kmer larger than 32nt
    @test_throws ArgumentError DNAMer(random_dna(33)) # no kmer larger than 32nt
    @test_throws ArgumentError RNAMer(
                      RNA_A, RNA_C, RNA_G, RNA_U, # no kmer larger than 32nt
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U)
    @test_throws ArgumentError DNAMer(
                      DNA_A, DNA_C, DNA_G, DNA_T, # no kmer larger than 32nt
                      DNA_A, DNA_C, DNA_G, DNA_T,
                      DNA_A, DNA_C, DNA_G, DNA_T,
                      DNA_A, DNA_C, DNA_G, DNA_T,
                      DNA_A, DNA_C, DNA_G, DNA_T,
                      DNA_A, DNA_C, DNA_G, DNA_T,
                      DNA_A, DNA_C, DNA_G, DNA_T,
                      DNA_A, DNA_C, DNA_G, DNA_T,
                      DNA_A, DNA_C, DNA_G, DNA_T)
    
    @test_throws ArgumentError BigDNAMer{0}(UInt128(0)) # 0-mers not allowed
    @test_throws ArgumentError BigRNAMer{0}(UInt128(0)) # 0-mers not allowed
    @test_throws ArgumentError BigRNAMer(RNA_A, RNA_C, RNA_G, RNA_N, RNA_U) # no Ns in kmers
    @test_throws ArgumentError BigDNAMer(DNA_A, DNA_C, DNA_G, DNA_N, DNA_T) # no Ns in kmers
    @test_throws ArgumentError BigRNAMer(rna"ACGNU")# no Ns in kmers
    @test_throws ArgumentError BigDNAMer(dna"ACGNT") # no Ns in kmers
    @test_throws MethodError BigDNAMer(RNA_A, DNA_A) # no mixing of RNA and DNA
    @test_throws ArgumentError BigRNAMer(random_rna(65)) # no kmer larger than 64nt
    @test_throws ArgumentError BigDNAMer(random_dna(65)) # no kmer larger than 64nt
    @test_throws ArgumentError BigRNAMer(
                      RNA_A, RNA_C, RNA_G, RNA_U, # no kmer larger than 64nt
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U,
                      RNA_A, RNA_C, RNA_G, RNA_U)
    @test_throws ArgumentError BigDNAMer(
                      DNA_A, DNA_C, DNA_G, DNA_T, # no kmer larger than 64nt
                      DNA_A, DNA_C, DNA_G, DNA_T,
                      DNA_A, DNA_C, DNA_G, DNA_T,
                      DNA_A, DNA_C, DNA_G, DNA_T,
                      DNA_A, DNA_C, DNA_G, DNA_T,
                      DNA_A, DNA_C, DNA_G, DNA_T,
                      DNA_A, DNA_C, DNA_G, DNA_T,
                      DNA_A, DNA_C, DNA_G, DNA_T,
                      DNA_A, DNA_C, DNA_G, DNA_T,
                      DNA_A, DNA_C, DNA_G, DNA_T,
                      DNA_A, DNA_C, DNA_G, DNA_T,
                      DNA_A, DNA_C, DNA_G, DNA_T,
                      DNA_A, DNA_C, DNA_G, DNA_T,
                      DNA_A, DNA_C, DNA_G, DNA_T,
                      DNA_A, DNA_C, DNA_G, DNA_T,
                      DNA_A, DNA_C, DNA_G, DNA_T,
                      DNA_A, DNA_C, DNA_G, DNA_T,
                      DNA_A, DNA_C, DNA_G, DNA_T)

    @testset "From strings" begin
        @test DNAMer("ACTG") == DNAMer(LongDNASeq("ACTG"))
        @test RNAMer("ACUG") == RNAMer(LongRNASeq("ACUG"))
        @test BigDNAMer("ACTG") == BigDNAMer(LongDNASeq("ACTG"))
        @test BigRNAMer("ACUG") == BigRNAMer(LongRNASeq("ACUG"))

        # N is not allowed in Kmers
        @test_throws Exception DNAMmer("ACGTNACGT")
        @test_throws Exception RNAMer("ACGUNACGU")
        @test_throws Exception BigDNAMmer("ACGTNACGT")
        @test_throws Exception BigRNAMer("ACGUNACGU")

        # Test string literals
        @test mer"ACTG"dna == DNAMer(LongDNASeq("ACTG"))
        @test isa(mer"ACGT"dna, DNAMer{4})
        @test_throws LoadError eval(:(mer"ACGN"dna))
        @test_throws LoadError eval(:(mer"ACG-"dna))
        
        @test bigmer"ACTG"dna == BigDNAMer(LongDNASeq("ACTG"))
        @test isa(bigmer"ACGT"dna, BigDNAMer{4})
        @test_throws LoadError eval(:(bigmer"ACGN"dna))
        @test_throws LoadError eval(:(bigmer"ACG-"dna))
    end
    
    @testset "Capacity" begin
        @test BioSequences.capacity(DNAMer(random_dna_kmer(10))) == 32
        @test BioSequences.capacity(RNAMer(random_rna_kmer(10))) == 32
        @test BioSequences.capacity(DNAMer(random_dna_kmer(32))) == 32
        @test BioSequences.capacity(RNAMer(random_rna_kmer(32))) == 32
        
        @test BioSequences.capacity(BigDNAMer(random_dna_kmer(10))) == 64
        @test BioSequences.capacity(BigRNAMer(random_rna_kmer(10))) == 64
        @test BioSequences.capacity(BigDNAMer(random_dna_kmer(32))) == 64
        @test BioSequences.capacity(BigRNAMer(random_rna_kmer(32))) == 64
    end
    
    @testset "N unused" begin
        @test BioSequences.n_unused(DNAMer(random_dna_kmer(10))) == 22
        @test BioSequences.n_unused(RNAMer(random_rna_kmer(10))) == 22
        @test BioSequences.n_unused(DNAMer(random_dna_kmer(32))) == 0
        @test BioSequences.n_unused(RNAMer(random_rna_kmer(32))) == 0
        
        @test BioSequences.n_unused(BigDNAMer(random_dna_kmer(40))) == 24
        @test BioSequences.n_unused(BigRNAMer(random_rna_kmer(40))) == 24
        @test BioSequences.n_unused(BigDNAMer(random_dna_kmer(64))) == 0
        @test BioSequences.n_unused(BigRNAMer(random_rna_kmer(64))) == 0
    end
end
