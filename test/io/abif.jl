@testset "ABIF Reader" begin
    function check_abif_parse(filename)
        stream = open(ABIF.Reader, filename)

        for record in stream end
        for (a,b) in collect(stream[1]) end
        for (a,b) in getindex(stream, ABIF.get_tags(stream)) end
        @test typeof(stream) == ABIF.Reader{IOStream}
        @test ABIF.get_tags(stream)[1].name == "AEPt"
        @test ABIF.get_tags(stream, "DATA")[1].name == "DATA"
        @test length(stream["DATA"]) == 12
        @test length(stream[1]) == 1

        @test typeof(getindex(stream, ABIF.get_tags(stream))) == Dict{String,Any}
        @test ABIF.tagelements(stream, "DATA") == 12
    end

    get_bio_fmt_specimens()
    path = joinpath(dirname(dirname(@__FILE__)), "BioFmtSpecimens", "ABI")
    for specimen in YAML.load_file(joinpath(path, "index.yml"))
        valid = get(specimen, "valid", true)
        filepath = joinpath(path, specimen["filename"])
        if valid
            check_abif_parse(filepath)
        else
            @test_throws Exception check_abif_parse(filepath)
        end
    end
end
