@testset "FASTQ" begin
    @test isa(FASTQ.Record("1", dna"AA", UInt8[10, 11]), FASTQ.Record)
    @test isa(FASTQ.Record("1", "desc.", dna"AA", UInt8[10, 11]), FASTQ.Record)
    @test_throws ArgumentError FASTQ.Record("1", dna"AA", UInt8[10])

    output = IOBuffer()
    writer = FASTQ.Writer(output)
    write(writer, FASTQ.Record("1", dna"AN", UInt8[11, 25]))
    write(writer, FASTQ.Record("2", "high quality", dna"TGA", UInt8[40, 41, 45]))
    flush(writer)
    @test String(take!(output)) == """
    @1
    AN
    +
    ,:
    @2 high quality
    TGA
    +
    IJN
    """

    output = IOBuffer()
    writer = FASTQ.Writer(output, quality_header=true)
    write(writer, FASTQ.Record("1", dna"AN", UInt8[11, 25]))
    write(writer, FASTQ.Record("2", "high quality", dna"TGA", UInt8[40, 41, 45]))
    flush(writer)
    @test String(take!(output)) == """
    @1
    AN
    +1
    ,:
    @2 high quality
    TGA
    +2 high quality
    IJN
    """

    @testset "Record" begin
        record = FASTQ.Record()
        @test !isfilled(record)

        record = FASTQ.Record("""
        @SRR1238088.1.1 HWI-ST499:111:D0G94ACXX:1:1101:1173:2105
        AAGCTCATGACCCGTCTTACCTACACCCTTGACGAGATCGAAGGA
        +SRR1238088.1.1 HWI-ST499:111:D0G94ACXX:1:1101:1173:2105
        @BCFFFDFHHHHHJJJIJIJJIJJJJJJJJIJJJJIIIJJJIJJJ
        """)
        @test isfilled(record)
        @test FASTQ.hasidentifier(record) == hasseqname(record) == true
        @test FASTQ.identifier(record) == seqname(record) == "SRR1238088.1.1"
        @test FASTQ.hasdescription(record)
        @test FASTQ.description(record) == "HWI-ST499:111:D0G94ACXX:1:1101:1173:2105"
        @test FASTQ.hassequence(record) == hassequence(record) == true
        @test FASTQ.sequence(DNASequence, record) == dna"AAGCTCATGACCCGTCTTACCTACACCCTTGACGAGATCGAAGGA"
        @test FASTQ.sequence(record) == sequence(record) == dna"AAGCTCATGACCCGTCTTACCTACACCCTTGACGAGATCGAAGGA"
        @test FASTQ.sequence(String, record) == "AAGCTCATGACCCGTCTTACCTACACCCTTGACGAGATCGAAGGA"
        @test FASTQ.hasquality(record)
        @test FASTQ.quality(record) == b"@BCFFFDFHHHHHJJJIJIJJIJJJJJJJJIJJJJIIIJJJIJJJ" .- 33

        record = FASTQ.Record("""
        @SRR1238088.1.1
        AAGCTCATGACCCGTCTTACCTACACCCTTGACGAGATCGAAGGA
        +
        @BCFFFDFHHHHHJJJIJIJJIJJJJJJJJIJJJJIIIJJJIJJJ
        """)
        @test isfilled(record)
        @test !FASTQ.hasdescription(record)
    end

    function test_records(rs1, rs2)
        if length(rs1) != length(rs2)
            return false
        end
        for (r1, r2) in zip(rs1, rs2)
            if FASTQ.identifier(r1) != FASTQ.identifier(r2) ||
               FASTQ.sequence(r1)   != FASTQ.sequence(r2)   ||
               FASTQ.quality(r1)    != FASTQ.quality(r2)
                return false
            end
        end
        return true
    end

    function test_fastq_parse(filename, valid)
        # Reading from a reader
        reader = open(FASTQ.Reader, filename)
        @test eltype(reader) == FASTQ.Record
        if valid
            for record in reader end
            @test true  # no error
            @test close(reader) === nothing
        else
            @test_throws Exception begin
                for record in reader end
            end
            return
        end

        # in-place parsing
        reader = open(FASTQ.Reader, filename)
        record = eltype(reader)()
        try
            while true
                read!(reader, record)
            end
        catch ex
            close(reader)
            if !isa(ex, EOFError)
                rethrow()
            end
        end

        # Check round trip
        output = IOBuffer()
        writer = FASTQ.Writer(output)
        expected_entries = FASTQ.Record[]
        for record in open(FASTQ.Reader, filename)
            write(writer, record)
            push!(expected_entries, record)
        end
        flush(writer)

        seekstart(output)
        read_entries = FASTQ.Record[]
        for record in FASTQ.Reader(output)
            push!(read_entries, record)
        end

        return test_records(expected_entries, read_entries)
    end

    get_bio_fmt_specimens()
    path = joinpath(dirname(@__FILE__), "..", "BioFmtSpecimens", "FASTQ")
    for specimen in YAML.load_file(joinpath(path, "index.yml"))
        tags = split(get(specimen, "tags", ""))
        valid = get(specimen, "valid", true)
        # currently unsupported features
        if any(t ∈ tags for t in ["gaps", "rna", "comments", "linewrap"])
            continue
        end
        filename = specimen["filename"]
        test_fastq_parse(joinpath(path, filename), valid)
    end

    @testset "invalid quality encoding" begin
        # Sanger full range (note escape characters before '$' and '\')
        record = FASTQ.Record("""
        @FAKE0001 Original version has PHRED scores from 0 to 93 inclusive (in that order)
        ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC
        +
        !"#\$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~
        """)

        # the range is not enough in these encodings
        for encoding in (:solexa, :illumina13, :illumina15)
            @test_throws ErrorException FASTQ.quality(record, encoding)
        end

        # the range is enough in these encodings
        for encoding in (:sanger, :illumina18)
            @test FASTQ.quality(record, encoding) == collect(0:93)
        end
    end

    @testset "fill ambiguous nucleotides" begin
        input = IOBuffer("""
        @seq1
        ACGTNRacgtnr
        +
        BBBB##AAAA##
        """)
        @test FASTQ.sequence(first(FASTQ.Reader(input, fill_ambiguous=nothing))) == dna"ACGTNRACGTNR"
        seekstart(input)
        @test FASTQ.sequence(first(FASTQ.Reader(input, fill_ambiguous=DNA_A)))   == dna"ACGTAAACGTAA"
        seekstart(input)
        @test FASTQ.sequence(first(FASTQ.Reader(input, fill_ambiguous=DNA_G)))   == dna"ACGTGGACGTGG"
        seekstart(input)
        @test FASTQ.sequence(first(FASTQ.Reader(input, fill_ambiguous=DNA_N)))   == dna"ACGTNNACGTNN"
        seekstart(input)
        @test FASTQ.sequence(GeneralSequence{DNAAlphabet{2}}, first(FASTQ.Reader(input, fill_ambiguous=DNA_A))) == dna"ACGTAAACGTAA"
    end
end

@testset "Quality scores" begin
    @testset "Decoding base quality scores" begin
        function test_decode(encoding, values, expected)
            result = Array{Int8}(undef, length(expected))
            FASTQ.decode_quality_string!(encoding, values, result, 1, length(result))
            @test result == expected
        end

        test_decode(FASTQ.SANGER_QUAL_ENCODING,
                    UInt8['!', '#', '$', '%', '&', 'I', '~'],
                    Int8[0, 2, 3, 4, 5, 40, 93])

        test_decode(FASTQ.SOLEXA_QUAL_ENCODING,
                    UInt8[';', 'B', 'C', 'D', 'E', 'h', '~'],
                    Int8[-5, 2, 3, 4, 5, 40, 62])

        test_decode(FASTQ.ILLUMINA13_QUAL_ENCODING,
                    UInt8['@', 'B', 'C', 'D', 'E', 'h', '~'],
                    Int8[0, 2, 3, 4, 5, 40, 62])

        test_decode(FASTQ.ILLUMINA15_QUAL_ENCODING,
                    UInt8['C', 'D', 'E', 'h', '~'],
                    Int8[3, 4, 5, 40, 62])

        test_decode(FASTQ.ILLUMINA18_QUAL_ENCODING,
                    UInt8['!', '#', '$', '%', '&', 'I', '~'],
                    Int8[0, 2, 3, 4, 5, 40, 93])
    end

    @testset "Encoding base quality scores" begin
        function test_encode(encoding, values, expected)
            result = Array{UInt8}(undef, length(expected))
            FASTQ.encode_quality_string!(encoding, values, result, 1, length(result))
            @test result == expected
        end

        test_encode(FASTQ.SANGER_QUAL_ENCODING,
                    Int8[0, 2, 3, 4, 5, 40, 93],
                    UInt8['!', '#', '$', '%', '&', 'I', '~'])

        test_encode(FASTQ.SOLEXA_QUAL_ENCODING,
                    Int8[-5, 2, 3, 4, 5, 40, 62],
                    UInt8[';', 'B', 'C', 'D', 'E', 'h', '~'])

        test_encode(FASTQ.ILLUMINA13_QUAL_ENCODING,
                    Int8[0, 2, 3, 4, 5, 40, 62],
                    UInt8['@', 'B', 'C', 'D', 'E', 'h', '~'])

        test_encode(FASTQ.ILLUMINA15_QUAL_ENCODING,
                    Int8[3, 4, 5, 40, 62],
                    UInt8['C', 'D', 'E', 'h', '~'])

        test_encode(FASTQ.ILLUMINA18_QUAL_ENCODING,
                    Int8[0, 2, 3, 4, 5, 40, 93],
                    UInt8['!', '#', '$', '%', '&', 'I', '~'])
    end
end
