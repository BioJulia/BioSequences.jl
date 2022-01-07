###
### Genetic Code
###
###
### Genetic code table and translator from RNA to amino acid sequence.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

const XNA = Union{DNA, RNA}
function unambiguous_codon(a::XNA, b::XNA, c::XNA)
    @inbounds begin
    bits = twobitnucs[reinterpret(UInt8, a) + 0x01] << 4 |
    twobitnucs[reinterpret(UInt8, b) + 0x01] << 2 |
    twobitnucs[reinterpret(UInt8, c) + 0x01]
    end
    #reinterpret(RNACodon, bits % UInt64)
    return bits % UInt64
end

# A genetic code is a table mapping RNA 3-mers (i.e. RNAKmer{3}) to AminoAcids.
"Type representing a Genetic Code"
struct GeneticCode <: AbstractDict{UInt64, AminoAcid}
    name::String
    tbl::NTuple{64, AminoAcid}
end

###
### Basic Functions
###

function Base.getindex(code::GeneticCode, codon::UInt64)
    return @inbounds code.tbl[codon + one(UInt64)]
end

Base.copy(code::GeneticCode) = GeneticCode(copy(code.name), copy(code.tbl))
Base.length(code::GeneticCode) = 64

Base.show(io::IO, code::GeneticCode) = print(io, code.name)

function Base.show(io::IO, ::MIME"text/plain", code::GeneticCode)
    print(io, code.name)
    rna = rna"ACGU"
    for x in rna, y in rna
        println(io)
        print(io, "  ")
        for z in rna
            codon = unambiguous_codon(x, y, z)
            aa = code[codon]
            print(io, x, y, z, ": ", aa)
            if z != RNA_U
                print(io, "    ")
            end
        end
    end
end

###
### Iterating through genetic code
###


function Base.iterate(code::GeneticCode, x = UInt64(0))
    if x > UInt64(0b111111)
        return nothing
    else
        return (x, @inbounds code[x]), x + 1
    end
end

###
### Default genetic codes
###

struct TransTables
    tables::Dict{Int,GeneticCode}
    bindings::Dict{Int,Symbol}
    function TransTables()
        return new(Dict(), Dict())
    end
end

Base.getindex(trans::TransTables, key::Integer) = trans.tables[Int(key)]

function Base.show(io::IO, trans::TransTables)
    print(io, "Translation Tables:")
    ids = sort(collect(keys(trans.tables)))
    for id in ids
        println(io)
        print(io, lpad(id, 3), ". ")
        show(io, trans.tables[id])
        if haskey(trans.bindings, id)
            print(io, " (", trans.bindings[id], ")")
        end
    end
end

"""
Genetic code list of NCBI.

The standard genetic code is `ncbi_trans_table[1]` and others can be shown by
`show(ncbi_trans_table)`.
For more details, consult the next link:
http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes.
"""
const ncbi_trans_table = TransTables()

macro register_ncbi_gencode(id, bind, tbl)
    quote
        gencode = parse_gencode($tbl)
        const $(esc(bind)) = gencode
        ncbi_trans_table.tables[$id] = gencode
        ncbi_trans_table.bindings[$id] = Symbol($(string(bind)))
    end
end

function parse_gencode(s)
    name, _, aas, _, base1, base2, base3 = split(chomp(s), '\n')
    name = split(name, ' ', limit = 2)[2]  # drop number
    codearr = fill(AA_X, 4^3)
    @assert length(aas) == 73
    for i in 10:73
        aa = AminoAcid(aas[i])
        b1 = DNA(base1[i])
        b2 = DNA(base2[i])
        b3 = DNA(base3[i])
        codon = unambiguous_codon(b1, b2, b3)
        codearr[codon + one(UInt64)] = aa
    end
    return GeneticCode(name, NTuple{64, AminoAcid}(codearr))
end

# Genetic codes translation tables are taken from the NCBI taxonomy database.

@register_ncbi_gencode 1 standard_genetic_code """
1. The Standard Code

  AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = ---M---------------M---------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""

@register_ncbi_gencode 2 vertebrate_mitochondrial_genetic_code """
2. The Vertebrate Mitochondrial Code

  AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG
Starts = --------------------------------MMMM---------------M------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""

@register_ncbi_gencode 3 yeast_mitochondrial_genetic_code """
3. The Yeast Mitochondrial Code

  AAs  = FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = ----------------------------------MM----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""

@register_ncbi_gencode 4 mold_mitochondrial_genetic_code """
4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code

  AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = --MM---------------M------------MMMM---------------M------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""

@register_ncbi_gencode 5 invertebrate_mitochondrial_genetic_code """
5. The Invertebrate Mitochondrial Code

  AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG
Starts = ---M----------------------------MMMM---------------M------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""

@register_ncbi_gencode 6 ciliate_nuclear_genetic_code """
6. The Ciliate, Dasycladacean and Hexamita Nuclear Code

  AAs  = FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = -----------------------------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""

@register_ncbi_gencode 9 echinoderm_mitochondrial_genetic_code """
9. The Echinoderm and Flatworm Mitochondrial Code

  AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG
Starts = -----------------------------------M---------------M------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""

@register_ncbi_gencode 10 euplotid_nuclear_genetic_code """
10. The Euplotid Nuclear Code

  AAs  = FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = -----------------------------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""

@register_ncbi_gencode 11 bacterial_plastid_genetic_code """
11. The Bacterial, Archaeal and Plant Plastid Code

  AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = ---M---------------M------------MMMM---------------M------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""

@register_ncbi_gencode 12 alternative_yeast_nuclear_genetic_code """
12. The Alternative Yeast Nuclear Code

  AAs  = FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = -------------------M---------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""

@register_ncbi_gencode 13 ascidian_mitochondrial_genetic_code """
13. The Ascidian Mitochondrial Code

  AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG
Starts = ---M------------------------------MM---------------M------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""

@register_ncbi_gencode 14 alternative_flatworm_mitochondrial_genetic_code """
14. The Alternative Flatworm Mitochondrial Code

  AAs  = FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG
Starts = -----------------------------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""

@register_ncbi_gencode 16 chlorophycean_mitochondrial_genetic_code """
16. Chlorophycean Mitochondrial Code

  AAs  = FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = -----------------------------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""

@register_ncbi_gencode 21 trematode_mitochondrial_genetic_code """
21. Trematode Mitochondrial Code

  AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG
Starts = -----------------------------------M---------------M------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""

@register_ncbi_gencode 22 scenedesmus_obliquus_mitochondrial_genetic_code """
22. Scenedesmus obliquus Mitochondrial Code

  AAs  = FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = -----------------------------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""

@register_ncbi_gencode 23 thraustochytrium_mitochondrial_genetic_code """
23. Thraustochytrium Mitochondrial Code

  AAs  = FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = --------------------------------M--M---------------M------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""

@register_ncbi_gencode 24 pterobrachia_mitochondrial_genetic_code """
24. Pterobranchia Mitochondrial Code

  AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG
Starts = ---M---------------M---------------M---------------M------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""

@register_ncbi_gencode 25 candidate_division_sr1_genetic_code """
25. Candidate Division SR1 and Gracilibacteria Code

  AAs  = FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = ---M-------------------------------M---------------M------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""

###
### Translation
###

"""
    translate(seq, code=standard_genetic_code, allow_ambiguous_codons=true, alternative_start=false)

Translate an `LongRNA` or a `LongDNA` to an `LongAA`.

Translation uses genetic code `code` to map codons to amino acids. See
`ncbi_trans_table` for available genetic codes.
If codons in the given sequence cannot determine a unique amino acid, they
will be translated to `AA_X` if `allow_ambiguous_codons` is `true` and otherwise
result in an error. For organisms that utilize alternative start codons, one
can set `alternative_start=true`, in which case the first codon will always be
converted to a methionine.
"""
function translate(ntseq::LongNuc{N};
    code::GeneticCode = standard_genetic_code,
    allow_ambiguous_codons::Bool = true,
    alternative_start::Bool = false
) where N
    len = div((length(ntseq) % UInt) * 11, 32)
    translate!(LongAA(undef, len), ntseq; code = code,
    allow_ambiguous_codons = allow_ambiguous_codons, alternative_start = alternative_start)
end

function translate!(aaseq::LongAA,
    ntseq::LongNuc{2};#Sequence{<:Union{DNAAlphabet{2}, RNAAlphabet{2}}};
    code::GeneticCode = standard_genetic_code,
    allow_ambiguous_codons::Bool = true,
    alternative_start::Bool = false
)
    n_aa, remainder = divrem(length(ntseq) % UInt, 3)
    iszero(remainder) || error("LongRNA length is not divisible by three. Cannot translate.")
    resize!(aaseq, n_aa)
    @inbounds for i in 1:n_aa
        a = ntseq[3i-2]
        b = ntseq[3i-1]
        c = ntseq[3i]
        codon = unambiguous_codon(a, b, c)
        aaseq[i] = code[codon]
    end
    alternative_start && !isempty(aaseq) && (@inbounds aaseq[1] = AA_M)
    aaseq
end

function translate!(aaseq::LongAA,
    ntseq::LongNuc{4};#Sequence{<:Union{DNAAlphabet{4}, RNAAlphabet{4}}};
    code::GeneticCode = standard_genetic_code,
    allow_ambiguous_codons::Bool = true,
    alternative_start::Bool = false
)
    n_aa, remainder = divrem(length(ntseq) % UInt, 3)
    iszero(remainder) || error("LongRNA length is not divisible by three. Cannot translate.")
    resize!(aaseq, n_aa)
    @inbounds for i in 1:n_aa
        a = reinterpret(RNA, ntseq[3i-2])
        b = reinterpret(RNA, ntseq[3i-1])
        c = reinterpret(RNA, ntseq[3i])
        if isambiguous(a) | isambiguous(b) | isambiguous(c)
            aa = try_translate_ambiguous_codon(code, a, b, c)
            if aa === nothing
                if allow_ambiguous_codons
                    aa = AA_X
                else
                    error("codon ", a, b, c, " cannot be unambiguously translated")
                end
            end
            aaseq[i] = aa
        else
            aaseq[i] = code[unambiguous_codon(a, b, c)]
        end
    end
    alternative_start && !isempty(aaseq) && (@inbounds aaseq[1] = AA_M)
    aaseq
end


function try_translate_ambiguous_codon(code::GeneticCode,
                                       x::RNA,
                                       y::RNA,
                                       z::RNA)
    @inbounds if !isambiguous(x) & !isambiguous(y)
        # try to translate a codon `(x, y, RNA_N)`
        aa_a = code[unambiguous_codon(x, y, RNA_A)]
        aa_c = code[unambiguous_codon(x, y, RNA_C)]
        aa_g = code[unambiguous_codon(x, y, RNA_G)]
        aa_u = code[unambiguous_codon(x, y, RNA_U)]
        if aa_a == aa_c == aa_g == aa_u
            return aa_a
        end
    end

    found::Union{AminoAcid, Nothing} = nothing
    for (codon, aa) in code
        # TODO: Make this more tidy - maybe reuse the decode method for DNAAlph{4}.
        a = reinterpret(RNA, 0x1 << ((codon >> 4) & 0x3))
        b = reinterpret(RNA, 0x1 << ((codon >> 2) & 0x3))
        c = reinterpret(RNA, 0x1 << (codon & 0x3))
        @inbounds if (iscompatible(x, a) & iscompatible(y, b) & iscompatible(z, c))
            if found === nothing
                found = aa
            elseif aa != found
                return nothing
            end
        end
    end
    return found
end