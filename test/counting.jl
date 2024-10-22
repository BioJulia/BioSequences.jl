# Collection of BioSequence, LongSequence w. differing alphabets, and SubSequence
bio_seqs = BioSequence[
    SimpleSeq("UGAUCUGUAGU"),
    SimpleSeq("")
]

long_nucs = [
    dna"",
    dna"TAGAGGAC",
    dna"WSCGCTYCG",
    dna"ATGC-CAACC",
    dna"TAGCGA--TCCGAA-TCGAAGC",
    dna"TCG-GTWYCYGYA-KKKAKCK--KNAWSSSTAGTCYKMNNACYWS",
    rna"",
    rna"UAGUGCUGAG",
    rna"WSUACGYKMNAC",
    rna"UGAUCGUAGUAGUCGUCGUGUC",
    rna"UBBSSDMHGVBRHA-YNKYDKAWSDYWC-VDCDUKCY"
]

strip_type(::Type{<:DNAAlphabet}) = DNAAlphabet
strip_type(::Type{<:RNAAlphabet}) = RNAAlphabet

for i in 1:lastindex(long_nucs)
    sq = long_nucs[i]
    if all(iscertain, sq)
        push!(long_nucs, LongSequence{strip_type(typeof(Alphabet(sq))){2}}(sq))
    end
end

long_aa = [
    aa"",
    aa"ATYCISL",
    aa"-PYCOSL",
    aa"LDJFDODIVLJK",
    aa"OIDHEOU-JNFDEDWQ",
    aa"OIDFENU-JDFDWDWQA",
    aa"WQKPC--LELIFEJ--",
]

sub_nucs = [
    view(dna"", 1:0),
    view(rna"UAGUC", 1:0),
    view(LongDNA{2}("TAGTCGTAGGCTGCGTGATGATGATAGTCGTGTCGTAGTA"), 2:38),
    view(dna"TAGTGCTATK-GAMTTWGTWSGTA-GCCCTAS-GAGTCTTA-W", 19:42),
    view(dna"TCG-GTWYCYGYA-KKKAKCK--KNAWSSSTAGTCYKMNNACYWS", 4:45)
]

for seq in long_nucs
    for span in [3:19, 7:31, 1:length(seq)]
        last(span) > length(seq) && continue
        push!(sub_nucs, view(seq, span))
    end
end

all_nucs = append!(copy(bio_seqs), long_nucs, sub_nucs)
all_seqs = append!(BioSequence[], all_nucs, long_aa)

@testset failfast=true "Counting" begin
    # GC content
    for i in all_nucs
        gcc = gc_content(i)
        @test isempty(i) ? isnan(gcc) : isapprox(gcc, sum(isGC, i; init=0) / length(i))
        @test count(isGC, i) == sum(isGC, i; init=0)
    end

    # Matches / mismatches
    for i in 1:length(all_seqs)-1, j in i+1:length(all_seqs)
        a = all_seqs[i]
        b = all_seqs[j]

        # We do this gymnastics because we want to throw depwarns when seeing
        # sequences of differing lengths, and we can't throw depwarns during test.
        (a, b) = length(a) <  length(b) ? (a, b) : (b, a)
        b = if length(b) == length(a)
            b
        elseif b isa SimpleSeq
            # Can't take view of a SimpleSeq
            continue
        else
            view(b, 1:length(a))
        end
        @test mismatches(a, b) == sum(splat(!=), zip(a, b); init=0)
        @test matches(a, b) == sum(splat(==), zip(a, b); init=0)
    end

    # Gaps
    for i in all_seqs
        @test count(isgap, i) == sum(isgap, i; init=0)
    end

    # Ambiguous
    for i in all_seqs
        @test count(isambiguous, i) == sum(isambiguous, i; init=0)
    end

    # Certain
    for i in all_seqs
        @test count(iscertain, i) == sum(iscertain, i; init=0)
    end

    # Count symbol
    for seq in all_seqs
        for symbol in Alphabet(seq)
            @test count(==(symbol), seq) == sum(i == symbol for i in seq; init=0)
        end
    end
    @test_throws EncodeError count(==(DNA_A), rna"UAG")
    @test_throws EncodeError count(==(DNA_M), LongDNA{2}(dna"TAG"))
    @test_throws EncodeError count(==(AA_C), dna"TAG")
    @test_throws EncodeError count(==(RNA_U), dna"TAG")
end
