@testset "FASTA" begin
    @testset "Record" begin
        record = FASTA.Record()
        @test !isfilled(record)
        @test_throws ArgumentError FASTA.identifier(record)

        record = FASTA.Record(">foo\nACGT\n")
        @test isfilled(record)
        @test hasseqname(record)
        @test FASTA.hasidentifier(record)
        @test seqname(record) == FASTA.identifier(record) == "foo"
        @test !FASTA.hasdescription(record)
        @test_throws BioCore.Exceptions.MissingFieldException FASTA.description(record)
        @test hassequence(record)
        @test FASTA.hassequence(record)
        @test FASTA.sequence(record) == dna"ACGT"
        @test FASTA.sequence(record, 2:3) == dna"CG"
        @test FASTA.sequence(String, record) == "ACGT"
        @test FASTA.sequence(String, record, 2:3) == "CG"
        @test record == FASTA.Record(">foo\nACGT\n")

        record = FASTA.Record("""
        >CYS1_DICDI fragment
        SCWSFSTTGNVEGQHFISQNKL
        VSLSEQNLVDCDHECMEYEGE
        """)
        @test isfilled(record)
        @test FASTA.identifier(record) == "CYS1_DICDI"
        @test FASTA.description(record) == "fragment"
        @test FASTA.sequence(record) == aa"SCWSFSTTGNVEGQHFISQNKLVSLSEQNLVDCDHECMEYEGE"
        @test FASTA.sequence(record, 10:15) == aa"NVEGQH"
    end

    output = IOBuffer()
    writer = FASTA.Writer(output, 5)
    write(writer, FASTA.Record("seq1", dna"TTA"))
    write(writer, FASTA.Record("seq2", "some description", dna"ACGTNN"))
    flush(writer)
    @test String(take!(output)) == """
    >seq1
    TTA
    >seq2 some description
    ACGTN
    N
    """

    reader = FASTA.Reader(IOBuffer(
    """
    >seqA some description
    QIKDLLVSSSTDLDTTLKMK
    ILELPFASGDLSM
    >seqB
    VLMALGMTDLFIPSANLTG*
    """))
    record = FASTA.Record()
    @test read!(reader, record) === record
    @test FASTA.identifier(record) == "seqA"
    @test FASTA.description(record) == "some description"
    @test FASTA.sequence(record) == aa"QIKDLLVSSSTDLDTTLKMKILELPFASGDLSM"
    @test read!(reader, record) === record
    @test FASTA.identifier(record) == "seqB"
    @test !FASTA.hasdescription(record)
    @test FASTA.sequence(record) == aa"VLMALGMTDLFIPSANLTG*"

    function test_fasta_parse(filename, valid)
        # Reading from a stream
        stream = open(FASTA.Reader, filename)
        @test eltype(stream) == FASTA.Record
        if valid
            for seqrec in stream end
            @test true  # no error
            @test close(stream) === nothing
        else
            @test_throws Exception begin
                for seqrec in stream end
            end
            return
        end

        # in-place parsing
        stream = open(FASTA.Reader, filename)
        entry = eltype(stream)()
        while !eof(stream)
            read!(stream, entry)
        end

        # Check round trip
        output = IOBuffer()
        writer = FASTA.Writer(output, width = 60)
        outputB = IOBuffer()
        writerB = FASTA.Writer(outputB, width = -1)
        expected_entries = Any[]
        for seqrec in open(FASTA.Reader, filename)
            write(writer, seqrec)
            write(writerB, seqrec)
            push!(expected_entries, seqrec)
        end
        flush(writer)
        flush(writerB)

        seekstart(output)
        seekstart(outputB)
        read_entries = FASTA.Record[]
        read_entriesB = FASTA.Record[]
        for seqrec in FASTA.Reader(output)
            push!(read_entries, seqrec)
        end
        for seqrec in FASTA.Reader(outputB)
            push!(read_entriesB, seqrec)
        end
        @test expected_entries == read_entries
        @test expected_entries == read_entriesB
    end

    get_bio_fmt_specimens()
    path = joinpath(dirname(@__FILE__), "..", "BioFmtSpecimens", "FASTA")
    for specimen in YAML.load_file(joinpath(path, "index.yml"))
        tags = specimen["tags"]
        valid = get(specimen, "valid", true)
        if occursin("comments", tags)
            # currently comments are not supported
            continue
        end
        test_fasta_parse(joinpath(path, specimen["filename"]), valid)
    end

    @testset "Faidx" begin
        fastastr = """
        >chr1
        CCACACCACACCCACACACC
        >chr2
        ATGCATGCATGCAT
        GCATGCATGCATGC
        >chr3
        AAATAGCCCTCATGTACGTCTCCTCCAAGCCCTGTTGTCTCTTACCCGGA
        TGTTCAACCAAAAGCTACTTACTACCTTTATTTTATGTTTACTTTTTATA
        >chr4
        TACTT
        """
        # generated with `samtools faidx`
        faistr = """
        chr1	20	6	20	21
        chr2	28	33	14	15
        chr3	100	69	50	51
        chr4	5	177	5	6
        """
        mktempdir() do dir
            filepath = joinpath(dir, "test.fa")
            write(filepath, fastastr)
            write(filepath * ".fai", faistr)
            open(FASTA.Reader, filepath, index=filepath * ".fai") do reader
                chr3 = reader["chr3"]
                @test FASTA.identifier(chr3) == "chr3"
                @test FASTA.sequence(chr3) == dna"""
                AAATAGCCCTCATGTACGTCTCCTCCAAGCCCTGTTGTCTCTTACCCGGA
                TGTTCAACCAAAAGCTACTTACTACCTTTATTTTATGTTTACTTTTTATA
                """

                chr2 = reader["chr2"]
                @test FASTA.identifier(chr2) == "chr2"
                @test FASTA.sequence(chr2) == dna"""
                ATGCATGCATGCAT
                GCATGCATGCATGC
                """

                chr4 = reader["chr4"]
                @test FASTA.identifier(chr4) == "chr4"
                @test FASTA.sequence(chr4) == dna"""
                TACTT
                """

                chr1 = reader["chr1"]
                @test FASTA.identifier(chr1) == "chr1"
                @test FASTA.sequence(chr1) == dna"""
                CCACACCACACCCACACACC
                """

                @test_throws ArgumentError reader["chr5"]
            end
        end

        # invalid index
        @test_throws ArgumentError FASTA.Reader(IOBuffer(fastastr), index=π)
    end

    @testset "append" begin
        intempdir() do
            filepath = "test.fa"
            writer = open(FASTA.Writer, filepath)
            write(writer, FASTA.Record("seq1", dna"AAA"))
            close(writer)
            writer = open(FASTA.Writer, filepath, append=true)
            write(writer, FASTA.Record("seq2", dna"CCC"))
            close(writer)
            seqs = open(collect, FASTA.Reader, filepath)
            @test length(seqs) == 2
        end
    end
end
