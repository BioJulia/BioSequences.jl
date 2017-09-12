var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#BioSequences.jl-1",
    "page": "Home",
    "title": "BioSequences.jl",
    "category": "section",
    "text": "_Biological sequences for julia_Latest release:(Image: Latest Release) (Image: BioSequences) (Image: BioSequences) (Image: License) (Image: ) (Image: BioJulia maintainer: bicycle1885) (Image: BioJulia maintainer: Ward9250)Development status:(Image: Build Status) (Image: Build status) (Image: codecov) (Image: )"
},

{
    "location": "index.html#Description-1",
    "page": "Home",
    "title": "Description",
    "category": "section",
    "text": "BioSequences.jl provides DNA, RNA and amino acid sequence data types for the julia language, with a comprehensive set of methods for common operations and IO of major sequence data formats.   "
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "Install BioSequences from the Julia REPL:julia> Pkg.add(\"BioSequences\")If you are interested in the cutting edge of the development, please check out the master branch to try new features before release."
},

{
    "location": "symbols.html#",
    "page": "Biological Symbols",
    "title": "Biological Symbols",
    "category": "page",
    "text": "CurrentModule = BioSequences\nDocTestSetup = quote\n    using BioSequences\nend"
},

{
    "location": "symbols.html#Biological-symbols-1",
    "page": "Biological Symbols",
    "title": "Biological symbols",
    "category": "section",
    "text": "The BioSequences module reexports the biological symbol (character) types that are provided by BioSymbols.jl:Type Meaning\nDNA DNA nucleotide\nRNA RNA nucleotide\nAminoAcid Amino acidThese symbols are elements of biological sequences, just as characters are elements of strings. See sections beginning from Introduction to the sequence data-types section for details."
},

{
    "location": "symbols.html#DNA-and-RNA-nucleotides-1",
    "page": "Biological Symbols",
    "title": "DNA and RNA nucleotides",
    "category": "section",
    "text": "Set of nucleotide symbols in BioSequences.jl covers IUPAC nucleotide base plus a gap symbol:Symbol Constant Meaning\n'A' DNA_A / RNA_A A; Adenine\n'C' DNA_C / RNA_C C; Cytosine\n'G' DNA_G / RNA_G G; Guanine\n'T' DNA_T T; Thymine (DNA only)\n'U' RNA_U U; Uracil (RNA only)\n'M' DNA_M / RNA_M A or C\n'R' DNA_R / RNA_R A or G\n'W' DNA_W / RNA_W A or T/U\n'S' DNA_S / RNA_S C or G\n'Y' DNA_Y / RNA_Y C or T/U\n'K' DNA_K / RNA_K G or T/U\n'V' DNA_V / RNA_V A or C or G; not T/U\n'H' DNA_H / RNA_H A or C or T; not G\n'D' DNA_D / RNA_D A or G or T/U; not C\n'B' DNA_B / RNA_B C or G or T/U; not A\n'N' DNA_N / RNA_N A or C or G or T/U\n'-' DNA_Gap / RNA_Gap Gap (none of the above)http://www.insdc.org/documents/feature_table.html#7.4.1Symbols are accessible as constants with DNA_ or RNA_ prefix:julia> DNA_A\nDNA_A\n\njulia> DNA_T\nDNA_T\n\njulia> RNA_U\nRNA_U\n\njulia> DNA_Gap\nDNA_Gap\n\njulia> typeof(DNA_A)\nBioSymbols.DNA\n\njulia> typeof(RNA_A)\nBioSymbols.RNA\nSymbols can be constructed by converting regular characters:julia> convert(DNA, 'C')\nDNA_C\n\njulia> convert(DNA, 'C') === DNA_C\ntrue\nEvery nucleotide is encoded using the lower 4 bits of a byte. An unambiguous nucleotide has only one set bit and the other bits are unset. The table below summarizes all unambiguous nucleotides and their corresponding bits. An ambiguous nucleotide is the bitwise OR of unambiguous nucleotides that the ambiguous nucleotide can take. For example, DNA_R (meaning the nucleotide is either DNA_A or DNA_G) is encoded as 0101 because 0101 is the bitwise OR of 0001 (DNA_A) and 0100 (DNA_G). The gap symbol is always 0000.NucleicAcid Bits\nDNA_A, RNA_A 0001\nDNA_C, RNA_C 0010\nDNA_G, RNA_G 0100\nDNA_T, RNA_U 1000The next examples demonstrate bit operations of DNA:julia> bits(reinterpret(UInt8, DNA_A))\n\"00000001\"\n\njulia> bits(reinterpret(UInt8, DNA_G))\n\"00000100\"\n\njulia> bits(reinterpret(UInt8, DNA_R))\n\"00000101\"\n\njulia> bits(reinterpret(UInt8, DNA_B))\n\"00001110\"\n\njulia> ~DNA_A\nDNA_B\n\njulia> DNA_A | DNA_G\nDNA_R\n\njulia> DNA_R & DNA_B\nDNA_G\n"
},

{
    "location": "symbols.html#Amino-acids-1",
    "page": "Biological Symbols",
    "title": "Amino acids",
    "category": "section",
    "text": "Set of amino acid symbols also covers IUPAC amino acid symbols plus a gap symbol:Symbol Constant Meaning\n'A' AA_A Alanine\n'R' AA_R Arginine\n'N' AA_N Asparagine\n'D' AA_D Aspartic acid (Aspartate)\n'C' AA_C Cysteine\n'Q' AA_Q Glutamine\n'E' AA_E Glutamic acid (Glutamate)\n'G' AA_G Glycine\n'H' AA_H Histidine\n'I' AA_I Isoleucine\n'L' AA_L Leucine\n'K' AA_K Lysine\n'M' AA_M Methionine\n'F' AA_F Phenylalanine\n'P' AA_P Proline\n'S' AA_S Serine\n'T' AA_T Threonine\n'W' AA_W Tryptophan\n'Y' AA_Y Tyrosine\n'V' AA_V Valine\n'O' AA_O Pyrrolysine\n'U' AA_U Selenocysteine\n'B' AA_B Aspartic acid or Asparagine\n'J' AA_J Leucine or Isoleucine\n'Z' AA_Z Glutamine or Glutamic acid\n'X' AA_X Any amino acid\n'*' AA_Term Termination codon\n'-' AA_Gap Gap (none of the above)http://www.insdc.org/documents/feature_table.html#7.4.3Symbols are accessible as constants with AA_ prefix:julia> AA_A\nAA_A\n\njulia> AA_Q\nAA_Q\n\njulia> AA_Term\nAA_Term\n\njulia> typeof(AA_A)\nBioSymbols.AminoAcid\nSymbols can be constructed by converting regular characters:julia> convert(AminoAcid, 'A')\nAA_A\n\njulia> convert(AminoAcid, 'P') === AA_P\ntrue\n"
},

{
    "location": "symbols.html#BioSymbols.alphabet",
    "page": "Biological Symbols",
    "title": "BioSymbols.alphabet",
    "category": "Function",
    "text": "alphabet(DNA)\n\nGet all symbols of DNA in sorted order.\n\nExamples\n\njulia> alphabet(DNA)\n(DNA_Gap, DNA_A, DNA_C, DNA_M, DNA_G, DNA_R, DNA_S, DNA_V, DNA_T, DNA_W, DNA_Y, DNA_H, DNA_K, DNA_D, DNA_B, DNA_N)\n\njulia> issorted(alphabet(DNA))\ntrue\n\n\n\n\nalphabet(RNA)\n\nGet all symbols of RNA in sorted order.\n\nExamples\n\njulia> alphabet(RNA)\n(RNA_Gap, RNA_A, RNA_C, RNA_M, RNA_G, RNA_R, RNA_S, RNA_V, RNA_U, RNA_W, RNA_Y, RNA_H, RNA_K, RNA_D, RNA_B, RNA_N)\n\njulia> issorted(alphabet(RNA))\ntrue\n\n\n\n\nalphabet(AminoAcid)\n\nGet all symbols of AminoAcid in sorted order.\n\nExamples\n\njulia> alphabet(AminoAcid)\n(AA_A, AA_R, AA_N, AA_D, AA_C, AA_Q, AA_E, AA_G, AA_H, AA_I, AA_L, AA_K, AA_M, AA_F, AA_P, AA_S, AA_T, AA_W, AA_Y, AA_V, AA_O, AA_U, AA_B, AA_J, AA_Z, AA_X, AA_Term, AA_Gap)\n\njulia> issorted(alphabet(AminoAcid))\ntrue\n\n\n\n\nGets the alphabet encoding of a given BioSequence.\n\n\n\n"
},

{
    "location": "symbols.html#BioSymbols.gap",
    "page": "Biological Symbols",
    "title": "BioSymbols.gap",
    "category": "Function",
    "text": "gap(DNA)\n\nReturn DNA_Gap.\n\n\n\ngap(RNA)\n\nReturn RNA_Gap.\n\n\n\ngap(AminoAcid)\n\nReturn AA_Gap.\n\n\n\n"
},

{
    "location": "symbols.html#BioSymbols.iscompatible",
    "page": "Biological Symbols",
    "title": "BioSymbols.iscompatible",
    "category": "Function",
    "text": "iscompatible(x::T, y::T) where T <: NucleicAcid\n\nTest if x and y are compatible with each other (i.e. x and y can be the same symbol).\n\nx and y must be the same type.\n\nExamples\n\njulia> iscompatible(DNA_A, DNA_A)\ntrue\n\njulia> iscompatible(DNA_C, DNA_N)  # DNA_N can be DNA_C\ntrue\n\njulia> iscompatible(DNA_C, DNA_R)  # DNA_R (A or G) cannot be DNA_C\nfalse\n\n\n\n\niscompatible(x::AminoAcid, y::AminoAcid)\n\nTest if x and y are compatible with each other.\n\nExamples\n\njulia> iscompatible(AA_A, AA_R)\nfalse\n\njulia> iscompatible(AA_A, AA_X)\ntrue\n\n\n\n\n"
},

{
    "location": "symbols.html#BioSymbols.isambiguous",
    "page": "Biological Symbols",
    "title": "BioSymbols.isambiguous",
    "category": "Function",
    "text": "isambiguous(nt::NucleicAcid)\n\nTest if nt is an ambiguous nucleotide.\n\n\n\nisambiguous(aa::AminoAcid)\n\nTest if aa is an ambiguous amino acid.\n\n\n\n"
},

{
    "location": "symbols.html#Other-functions-1",
    "page": "Biological Symbols",
    "title": "Other functions",
    "category": "section",
    "text": "alphabet\ngap\niscompatible\nisambiguous"
},

{
    "location": "sequences/sequences.html#",
    "page": "Overview",
    "title": "Overview",
    "category": "page",
    "text": ""
},

{
    "location": "sequences/sequences.html#Biological-sequences-1",
    "page": "Overview",
    "title": "Biological sequences",
    "category": "section",
    "text": "The BioSequences module provides representations and tools for manipulating nucleotide and amino acid sequences."
},

{
    "location": "sequences/sequences.html#Introduction-to-the-sequence-data-types-1",
    "page": "Overview",
    "title": "Introduction to the sequence data-types",
    "category": "section",
    "text": "Sequences in BioSequences.jl are more strictly typed than in many other libraries; elements in a sequence are typed as biological symbol instead of character or byte. They are special purpose types rather than simply strings and hence offer additional functionality that naive string types don't have. Though this strictness sacrifices some convenience, it also means you can always rely on a DNA sequence type to store DNA and nothing but DNA, without having to check, or deal with lowercase versus uppercase and so on. Strict separation of sequence types also means we are free to choose the most efficient representation. DNA and RNA sequences are encoded using either four bits per base (which is the default), or two bits per base. This makes them memory efficient and allows us to speed up many common operations and transformations, like nucleotide composition, reverse complement, and k-mer enumeration.The BioSequences provides three different sequence types: BioSequence, Kmer and ReferenceSequence. Each of these types is a subtype of an abstract type called Sequence and supports various string-like operations such as random access and iteration. Different sequence types have different features. In most situations, BioSequence type will do and is used as the default representation. But sometimes other types are much more preferable in terms of memory efficiency and computation performance.  Here is the summary table of these three types:Type Description Element type Mutability Allocation\nBioSequence{A<:Alphabet} general-purpose biological sequences DNA, RNA, Amino acids mutable heap\nKmer{T<:NucleicAcid,k} specialized for short nucleotide sequences DNA, RNA immutable stack / register\nReferenceSequence specialized for long reference genomes DNA immutable heapDetails of these different representations are explained in the following sections:BioSequence: General-purpose sequences\nKmer: Nucleic acid k-mers\nReferenceSequence: Reference sequences"
},

{
    "location": "sequences/bioseq.html#",
    "page": "BioSequence",
    "title": "BioSequence",
    "category": "page",
    "text": "CurrentModule = BioSequences\nDocTestSetup = quote\n    using BioSequences\nend"
},

{
    "location": "sequences/bioseq.html#General-purpose-sequences-1",
    "page": "BioSequence",
    "title": "General-purpose sequences",
    "category": "section",
    "text": "BioSequence{A} is a generic sequence type parameterized by an alphabet type A that defines the domain (or set) of biological symbols, and each alphabet has an associated symbol type. For example, AminoAcidAlphabet is associated with AminoAcid and hence an object of the BioSequence{AminoAcidAlphabet} type represents a sequence of amino acids.  Symbols from multiple alphabets can't be intermixed in one sequence type.The following table summarizes common sequence types that are defined in the BioSequences module:Type Symbol type Type alias\nBioSequence{DNAAlphabet{4}} DNA DNASequence\nBioSequence{RNAAlphabet{4}} RNA RNASequence\nBioSequence{AminoAcidAlphabet} AminoAcid AminoAcidSequence\nBioSequence{CharAlphabet} Char CharSequenceParameterized definition of the BioSequence{A} type is for the purpose of unifying the data structure and operations of any symbol type. In most cases, users don't have to care about it and can use type aliases listed above. However, the alphabet type fixes the internal memory encoding and plays an important role when optimizing performance of a program (see Using a more compact sequence representation section for low-memory encodings).  It also enables a user to define their own alphabet only by defining few numbers of methods. This is described in Defining a new alphabet section."
},

{
    "location": "sequences/bioseq.html#Constructing-sequences-1",
    "page": "BioSequence",
    "title": "Constructing sequences",
    "category": "section",
    "text": ""
},

{
    "location": "sequences/bioseq.html#Using-string-literals-1",
    "page": "BioSequence",
    "title": "Using string literals",
    "category": "section",
    "text": "Most immediately, sequence literals can be constructed using the string macros dna, rna, aa, and char:julia> dna\"TACGTANNATC\"\n11nt DNA Sequence:\nTACGTANNATC\n\njulia> rna\"AUUUGNCCANU\"\n11nt RNA Sequence:\nAUUUGNCCANU\n\njulia> aa\"ARNDCQEGHILKMFPSTWYVX\"\n21aa Amino Acid Sequence:\nARNDCQEGHILKMFPSTWYVX\n\njulia> char\"αβγδϵ\"\n5char Char Sequence:\nαβγδϵ\nHowever it should be noted that by default these sequence literals allocate the BioSequence object before the code containing the sequence literal is run. This means there may be occasions where your program does not behave as you first expect. For example consider the following code:julia> function foo()\n           s = dna\"CTT\"\n           push!(s, DNA_A)\n       end\nfoo (generic function with 1 method)\nDocTestSetup = quote\n    using BioSequences\n    function foo()\n        s = dna\"CTT\"d\n        push!(s, DNA_A)\n    end\nendYou might expect that every time you call foo, that a DNA sequence CTTA would be returned. You might expect that this is because every time foo is called, a new DNA sequence variable CTT is created, and and A nucleotide is pushed to it, and the result, CTTA is returned. In other words you might expect the following output:julia> foo()\n4nt DNA Sequence:\nCTTA\n\njulia> foo()\n4nt DNA Sequence:\nCTTA\n\njulia> foo()\n4nt DNA Sequence:\nCTTA\nHowever, this is not what happens, instead the following happens:DocTestSetup = quote\n    using BioSequences\n    function foo()\n        s = dna\"CTT\"s\n        push!(s, DNA_A)\n    end\nendjulia> foo()\n4nt DNA Sequence:\nCTTA\n\njulia> foo()\n5nt DNA Sequence:\nCTTAA\n\njulia> foo()\n6nt DNA Sequence:\nCTTAAA\nThe reason for this is because the sequence literal is allocated only once before the first time the function foo is called and run. Therefore, s in foo is always a reference to that one sequence that was allocated. So one sequence is created before foo is called, and then it is pushed to every time foo is called. Thus, that one allocated sequence grows with every call of foo.If you wanted foo to create a new sequence each time it is called, then you can add a flag to the end of the sequence literal to dictate behaviour: A flag of 's' means 'static': the sequence will be allocated before code is run, as is the default behaviour described above. However providing 'd' flag changes the behaviour: 'd' means 'dynamic': the sequence will be allocated at whilst the code is running, and not before. So to change foo so as it creates a new sequence each time it is called, simply add the 'd' flag to the sequence literal:DocTestSetup = quote\n    using BioSequences\nendjulia> function foo()\n           s = dna\"CTT\"d     # 'd' flag appended to the string literal.\n           push!(s, DNA_A)\n       end\nfoo (generic function with 1 method)\nNow every time foo is called, a new sequence CTT is created, and an A nucleotide is pushed to it:DocTestSetup = quote\n    using BioSequences\n    function foo()\n        s = dna\"CTT\"d\n        push!(s, DNA_A)\n    end\nendjulia> foo()\n4nt DNA Sequence:\nCTTA\n\njulia> foo()\n4nt DNA Sequence:\nCTTA\n\njulia> foo()\n4nt DNA Sequence:\nCTTA\nSo the take home message of sequence literals is this:Be careful when you are using sequence literals inside of functions, and inside the bodies of things like for loops. And if you use them and are unsure, use the  's' and 'd' flags to ensure the behaviour you get is the behaviour you intend."
},

{
    "location": "sequences/bioseq.html#Other-constructors-and-conversion-1",
    "page": "BioSequence",
    "title": "Other constructors and conversion",
    "category": "section",
    "text": "Sequences can also be constructed from strings or arrays of nucleotide or amino acid symbols using constructors or the convert function:julia> DNASequence(\"TTANC\")\n5nt DNA Sequence:\nTTANC\n\njulia> DNASequence([DNA_T, DNA_T, DNA_A, DNA_N, DNA_C])\n5nt DNA Sequence:\nTTANC\n\njulia> convert(DNASequence, [DNA_T, DNA_T, DNA_A, DNA_N, DNA_C])\n5nt DNA Sequence:\nTTANC\nUsing convert, these operations are reversible: sequences can be converted to strings or arrays:julia> convert(String, dna\"TTANGTA\")\n\"TTANGTA\"\n\njulia> convert(Vector{DNA}, dna\"TTANGTA\")\n7-element Array{BioSymbols.DNA,1}:\n DNA_T\n DNA_T\n DNA_A\n DNA_N\n DNA_G\n DNA_T\n DNA_A\nSequences can also be concatenated into longer sequences:julia> DNASequence(dna\"ACGT\", dna\"NNNN\", dna\"TGCA\")\n12nt DNA Sequence:\nACGTNNNNTGCA\n\njulia> dna\"ACGT\" * dna\"TGCA\"\n8nt DNA Sequence:\nACGTTGCA\n\njulia> repeat(dna\"TA\", 10)\n20nt DNA Sequence:\nTATATATATATATATATATA\n\njulia> dna\"TA\" ^ 10\n20nt DNA Sequence:\nTATATATATATATATATATA\nDespite being separate types, DNASequence and RNASequence can freely be converted between efficiently without copying the underlying data:julia> dna = dna\"TTANGTAGACCG\"\n12nt DNA Sequence:\nTTANGTAGACCG\n\njulia> rna = convert(RNASequence, dna)\n12nt RNA Sequence:\nUUANGUAGACCG\n\njulia> dna.data === rna.data  # underlying data are same\ntrue\nA random sequence can be obtained by the randdnaseq, randrnaseq and randaaseq functions, which generate DNASequence, RNASequence and AminoAcidSequence, respectively. Generated sequences are composed of the standard symbols without ambiguity and gap. For example, randdnaseq(6) may generate dna\"TCATAG\" but never generates dna\"TNANAG\" or dna\"T-ATAG\".A translatable RNASequence can also be converted to an AminoAcidSequence using the translate function."
},

{
    "location": "sequences/bioseq.html#Indexing,-modifying-and-transformations-1",
    "page": "BioSequence",
    "title": "Indexing, modifying and transformations",
    "category": "section",
    "text": ""
},

{
    "location": "sequences/bioseq.html#Getindex-1",
    "page": "BioSequence",
    "title": "Getindex",
    "category": "section",
    "text": "Sequences for the most part behave like other vector or string types. They can be indexed using integers or ranges:julia> seq = dna\"ACGTTTANAGTNNAGTACC\"\n19nt DNA Sequence:\nACGTTTANAGTNNAGTACC\n\njulia> seq[5]\nDNA_T\n\njulia> seq[6:end]\n14nt DNA Sequence:\nTANAGTNNAGTACC\nNote that, indexing a biological sequence by range creates a subsequence of the original sequence. Unlike Arrays in the standard library, creating a subsequence is copy-free: a subsequence simply points to the original sequence data with its range. You may think that this is unsafe because modifying subsequences propagates to the original sequence, but this doesn't happen actually:julia> seq = dna\"AAAA\"    # create a sequence\n4nt DNA Sequence:\nAAAA\n\njulia> subseq = seq[1:2]  # create a subsequence from `seq`\n2nt DNA Sequence:\nAA\n\njulia> subseq[2] = DNA_T  # modify the second element of it\nDNA_T\n\njulia> subseq             # the subsequence is modified\n2nt DNA Sequence:\nAT\n\njulia> seq                # but the original sequence is not\n4nt DNA Sequence:\nAAAA\nThis is because modifying a sequence checks whether its underlying data are shared with other sequences under the hood. If and only if the data are shared, the subsequence creates a copy of itself. Any modifying operation does this check. This is called copy-on-write strategy and users don't need to care about it because it is transparent: If the user modifies a sequence with or subsequence, the job of managing and protecting the underlying data of sequences is handled for them."
},

{
    "location": "sequences/bioseq.html#Base.push!",
    "page": "BioSequence",
    "title": "Base.push!",
    "category": "Function",
    "text": "push!(collection, items...) -> collection\n\nInsert one or more items at the end of collection.\n\njulia> push!([1, 2, 3], 4, 5, 6)\n6-element Array{Int64,1}:\n 1\n 2\n 3\n 4\n 5\n 6\n\nUse append! to add all the elements of another collection to collection. The result of the preceding example is equivalent to append!([1, 2, 3], [4, 5, 6]).\n\n\n\n"
},

{
    "location": "sequences/bioseq.html#Base.pop!",
    "page": "BioSequence",
    "title": "Base.pop!",
    "category": "Function",
    "text": "pop!(collection, key[, default])\n\nDelete and return the mapping for key if it exists in collection, otherwise return default, or throw an error if default is not specified.\n\njulia> d = Dict(\"a\"=>1, \"b\"=>2, \"c\"=>3);\n\njulia> pop!(d, \"a\")\n1\n\njulia> pop!(d, \"d\")\nERROR: KeyError: key \"d\" not found\nStacktrace:\n [1] pop!(::Dict{String,Int64}, ::String) at ./dict.jl:539\n\njulia> pop!(d, \"e\", 4)\n4\n\n\n\npop!(collection) -> item\n\nRemove the last item in collection and return it.\n\njulia> A=[1, 2, 3, 4, 5, 6]\n6-element Array{Int64,1}:\n 1\n 2\n 3\n 4\n 5\n 6\n\njulia> pop!(A)\n6\n\njulia> A\n5-element Array{Int64,1}:\n 1\n 2\n 3\n 4\n 5\n\n\n\npop!(seq::BioSequence)\n\nRemove the symbol from the end of a biological sequence seq and return it. Returns a variable of eltype(seq).\n\n\n\n"
},

{
    "location": "sequences/bioseq.html#Base.shift!",
    "page": "BioSequence",
    "title": "Base.shift!",
    "category": "Function",
    "text": "shift!(collection) -> item\n\nRemove the first item from collection.\n\njulia> A = [1, 2, 3, 4, 5, 6]\n6-element Array{Int64,1}:\n 1\n 2\n 3\n 4\n 5\n 6\n\njulia> shift!(A)\n1\n\njulia> A\n5-element Array{Int64,1}:\n 2\n 3\n 4\n 5\n 6\n\n\n\n"
},

{
    "location": "sequences/bioseq.html#Base.unshift!",
    "page": "BioSequence",
    "title": "Base.unshift!",
    "category": "Function",
    "text": "unshift!(collection, items...) -> collection\n\nInsert one or more items at the beginning of collection.\n\njulia> unshift!([1, 2, 3, 4], 5, 6)\n6-element Array{Int64,1}:\n 5\n 6\n 1\n 2\n 3\n 4\n\n\n\nunshift!(seq, x)\n\nInsert a biological symbol x at the beginning of a biological sequence seq.\n\n\n\n"
},

{
    "location": "sequences/bioseq.html#Base.insert!",
    "page": "BioSequence",
    "title": "Base.insert!",
    "category": "Function",
    "text": "insert!(a::Vector, index::Integer, item)\n\nInsert an item into a at the given index. index is the index of item in the resulting a.\n\njulia> insert!([6, 5, 4, 2, 1], 4, 3)\n6-element Array{Int64,1}:\n 6\n 5\n 4\n 3\n 2\n 1\n\n\n\ninsert!(seq, i, x)\n\nInsert a biological symbol x into a biological sequence seq, at the given index i.\n\n\n\n"
},

{
    "location": "sequences/bioseq.html#Base.deleteat!-Tuple{BioSequences.BioSequence,Integer}",
    "page": "BioSequence",
    "title": "Base.deleteat!",
    "category": "Method",
    "text": "deleteat!(seq::BioSequence, i::Integer)\n\nDelete a biological symbol at a single position i in a biological sequence seq.\n\nModifies the input sequence.\n\n\n\n"
},

{
    "location": "sequences/bioseq.html#Base.append!",
    "page": "BioSequence",
    "title": "Base.append!",
    "category": "Function",
    "text": "append!(collection, collection2) -> collection.\n\nAdd the elements of collection2 to the end of collection.\n\njulia> append!([1],[2,3])\n3-element Array{Int64,1}:\n 1\n 2\n 3\n\njulia> append!([1, 2, 3], [4, 5, 6])\n6-element Array{Int64,1}:\n 1\n 2\n 3\n 4\n 5\n 6\n\nUse push! to add individual items to collection which are not already themselves in another collection. The result is of the preceding example is equivalent to push!([1, 2, 3], 4, 5, 6).\n\n\n\n"
},

{
    "location": "sequences/bioseq.html#Base.resize!",
    "page": "BioSequence",
    "title": "Base.resize!",
    "category": "Function",
    "text": "resize!(a::Vector, n::Integer) -> Vector\n\nResize a to contain n elements. If n is smaller than the current collection length, the first n elements will be retained. If n is larger, the new elements are not guaranteed to be initialized.\n\njulia> resize!([6, 5, 4, 3, 2, 1], 3)\n3-element Array{Int64,1}:\n 6\n 5\n 4\n\njulia> resize!([6, 5, 4, 3, 2, 1], 8)\n8-element Array{Int64,1}:\n 6\n 5\n 4\n 3\n 2\n 1\n 0\n 0\n\n\n\nresize!(seq, size)\n\nResize a biological sequence seq, to a given size.\n\n\n\n"
},

{
    "location": "sequences/bioseq.html#Setindex-and-modifying-DNA-sequences-1",
    "page": "BioSequence",
    "title": "Setindex and modifying DNA sequences",
    "category": "section",
    "text": "The biological symbol at a given locus in a biological sequence can be set using setindex:julia> seq = dna\"ACGTTTANAGTNNAGTACC\"\n19nt DNA Sequence:\nACGTTTANAGTNNAGTACC\n\njulia> seq[5] = DNA_A\nDNA_A\nIn addition, many other modifying operations are possible for biological sequences such as push!, pop!, and insert!, which should be familiar to people used to editing arrays.push!\npop!\nshift!\nunshift!\ninsert!\ndeleteat!(::BioSequences.BioSequence, ::Integer)\nappend!\nresize!Here are some examples:julia> seq = dna\"ACG\"\n3nt DNA Sequence:\nACG\n\njulia> push!(seq, DNA_T)\n4nt DNA Sequence:\nACGT\n\njulia> append!(seq, dna\"AT\")\n6nt DNA Sequence:\nACGTAT\n\njulia> deleteat!(seq, 2)\n5nt DNA Sequence:\nAGTAT\n\njulia> deleteat!(seq, 2:3)\n3nt DNA Sequence:\nAAT\n"
},

{
    "location": "sequences/bioseq.html#Base.reverse!-Tuple{BioSequences.BioSequence}",
    "page": "BioSequence",
    "title": "Base.reverse!",
    "category": "Method",
    "text": "reverse!(seq)\n\nReverse a biological sequence seq in place.\n\n\n\n"
},

{
    "location": "sequences/bioseq.html#BioSequences.complement!",
    "page": "BioSequence",
    "title": "BioSequences.complement!",
    "category": "Function",
    "text": "complement!(seq)\n\nMake a complement sequence of seq in place.\n\n\n\ncomplement!(seq)\n\nTransform seq into it's complement.\n\n\n\n"
},

{
    "location": "sequences/bioseq.html#BioSequences.reverse_complement!",
    "page": "BioSequence",
    "title": "BioSequences.reverse_complement!",
    "category": "Function",
    "text": "reverse_complement!(seq)\n\nMake a reversed complement sequence of seq in place.\n\nAmbiguous nucleotides are left as-is.\n\n\n\n"
},

{
    "location": "sequences/bioseq.html#BioSequences.ungap!",
    "page": "BioSequence",
    "title": "BioSequences.ungap!",
    "category": "Function",
    "text": "Remove gap characters from a sequence. Modifies the input sequence.\n\n\n\n"
},

{
    "location": "sequences/bioseq.html#Base.empty!",
    "page": "BioSequence",
    "title": "Base.empty!",
    "category": "Function",
    "text": "empty!(collection) -> collection\n\nRemove all elements from a collection.\n\njulia> A = Dict(\"a\" => 1, \"b\" => 2)\nDict{String,Int64} with 2 entries:\n  \"b\" => 2\n  \"a\" => 1\n\njulia> empty!(A);\n\njulia> A\nDict{String,Int64} with 0 entries\n\n\n\nempty!(seq)\n\nCompletely empty a biological sequence seq of nucleotides.\n\n\n\n"
},

{
    "location": "sequences/bioseq.html#Additional-transformations-1",
    "page": "BioSequence",
    "title": "Additional transformations",
    "category": "section",
    "text": "In addition to these basic modifying functions, other sequence transformations which are common in bioinformatics are also provided.reverse!(::BioSequences.BioSequence)\ncomplement!\nreverse_complement!\nungap!\nempty!Some examples:julia> seq = dna\"ACGTAT\"\n6nt DNA Sequence:\nACGTAT\n\njulia> reverse!(seq)\n6nt DNA Sequence:\nTATGCA\n\njulia> complement!(seq)\n6nt DNA Sequence:\nATACGT\n\njulia> reverse_complement!(seq)\n6nt DNA Sequence:\nACGTAT\nMany of these methods also have a version which makes a copy of the input sequence, so you get a modified copy, and don't alter the original sequence. Such methods are named the same, but without the exclamation mark. E.g. reverse instead of reverse!, and ungap instead of ungap!.  "
},

{
    "location": "sequences/bioseq.html#BioSequences.translate",
    "page": "BioSequence",
    "title": "BioSequences.translate",
    "category": "Function",
    "text": "translate(rna_seq, code=standard_genetic_code, allow_ambiguous_codons=true)\n\nTranslate an RNASequence to an AminoAcidSequence.\n\nTranslation uses genetic code code to map codons to amino acids. See ncbi_trans_table for available genetic codes. If codons in the given RNA sequence cannot determine a unique amino acid, they will be translated to AA_X if allow_ambiguous_codons is true and otherwise result in an error.\n\n\n\n"
},

{
    "location": "sequences/bioseq.html#BioSequences.ncbi_trans_table",
    "page": "BioSequence",
    "title": "BioSequences.ncbi_trans_table",
    "category": "Constant",
    "text": "Genetic code list of NCBI.\n\nThe standard genetic code is ncbi_trans_table[1] and others can be shown by show(ncbi_trans_table). For more details, consult the next link: http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes.\n\n\n\n"
},

{
    "location": "sequences/bioseq.html#Translation-1",
    "page": "BioSequence",
    "title": "Translation",
    "category": "section",
    "text": "Translation is a slightly more complex transformation for RNA Sequences and so we describe it here in more detail.The translate funtion translates a sequence of codons in a RNA sequence to a amino acid sequence besed on a genetic code mapping. The BioSequences module contains all NCBI defined genetic codes and they are registered in ncbi_trans_table.translate\nncbi_trans_tablejulia> ncbi_trans_table\nTranslation Tables:\n  1. The Standard Code (standard_genetic_code)\n  2. The Vertebrate Mitochondrial Code (vertebrate_mitochondrial_genetic_code)\n  3. The Yeast Mitochondrial Code (yeast_mitochondrial_genetic_code)\n  4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code (mold_mitochondrial_genetic_code)\n  5. The Invertebrate Mitochondrial Code (invertebrate_mitochondrial_genetic_code)\n  6. The Ciliate, Dasycladacean and Hexamita Nuclear Code (ciliate_nuclear_genetic_code)\n  9. The Echinoderm and Flatworm Mitochondrial Code (echinoderm_mitochondrial_genetic_code)\n 10. The Euplotid Nuclear Code (euplotid_nuclear_genetic_code)\n 11. The Bacterial, Archaeal and Plant Plastid Code (bacterial_plastid_genetic_code)\n 12. The Alternative Yeast Nuclear Code (alternative_yeast_nuclear_genetic_code)\n 13. The Ascidian Mitochondrial Code (ascidian_mitochondrial_genetic_code)\n 14. The Alternative Flatworm Mitochondrial Code (alternative_flatworm_mitochondrial_genetic_code)\n 16. Chlorophycean Mitochondrial Code (chlorophycean_mitochondrial_genetic_code)\n 21. Trematode Mitochondrial Code (trematode_mitochondrial_genetic_code)\n 22. Scenedesmus obliquus Mitochondrial Code (scenedesmus_obliquus_mitochondrial_genetic_code)\n 23. Thraustochytrium Mitochondrial Code (thraustochytrium_mitochondrial_genetic_code)\n 24. Pterobranchia Mitochondrial Code (pterobrachia_mitochondrial_genetic_code)\n 25. Candidate Division SR1 and Gracilibacteria Code (candidate_division_sr1_genetic_code)\nhttp://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes"
},

{
    "location": "sequences/bioseq.html#Site-counting-1",
    "page": "BioSequence",
    "title": "Site counting",
    "category": "section",
    "text": "BioSequences extends the Base.count method to provide some useful utilities for counting the number of sites in biological sequences."
},

{
    "location": "sequences/bioseq.html#BioSequences.Certain",
    "page": "BioSequence",
    "title": "BioSequences.Certain",
    "category": "Type",
    "text": "A Certain site describes a site where both of two aligned sites are not an ambiguity symbol or a gap.\n\n\n\n"
},

{
    "location": "sequences/bioseq.html#BioSequences.Gap",
    "page": "BioSequence",
    "title": "BioSequences.Gap",
    "category": "Type",
    "text": "An Gap site describes a site where either of two aligned sites are a gap symbol '-'.\n\n\n\n"
},

{
    "location": "sequences/bioseq.html#BioSequences.Ambiguous",
    "page": "BioSequence",
    "title": "BioSequences.Ambiguous",
    "category": "Type",
    "text": "An Ambiguous site describes a site where either of two aligned sites are an ambiguity symbol.\n\n\n\n"
},

{
    "location": "sequences/bioseq.html#BioSequences.Match",
    "page": "BioSequence",
    "title": "BioSequences.Match",
    "category": "Type",
    "text": "A Match site describes a site where two aligned nucleotides are the same biological symbol.\n\n\n\n"
},

{
    "location": "sequences/bioseq.html#BioSequences.Mismatch",
    "page": "BioSequence",
    "title": "BioSequences.Mismatch",
    "category": "Type",
    "text": "A Mismatch site describes a site where two aligned nucleotides are not the same biological symbol.\n\n\n\n"
},

{
    "location": "sequences/bioseq.html#Site-types-1",
    "page": "BioSequence",
    "title": "Site types",
    "category": "section",
    "text": "Different types of site can be counted. Each of these types is a concrete subtype of the abstract type Site:Certain\nGap\nAmbiguous\nMatch\nMismatch"
},

{
    "location": "sequences/bioseq.html#Base.count-methods-1",
    "page": "BioSequence",
    "title": "Base.count methods",
    "category": "section",
    "text": "The count method can be used with two sequences and a concrete subtype of Site:julia> count(Match, dna\"ATCGATCG\", dna\"AAGGTTCG\")\n5By providing a window and step size, counting can be done from within a sliding window:julia> count(Match, dna\"ATCGATCG\", dna\"AAGGTTCG\", 3, 1)\n6-element Array{IntervalTrees.IntervalValue{Int64,Int64},1}:\n IntervalTrees.IntervalValue{Int64,Int64}\n(1,3) => 1\n IntervalTrees.IntervalValue{Int64,Int64}\n(2,4) => 1\n IntervalTrees.IntervalValue{Int64,Int64}\n(3,5) => 1\n IntervalTrees.IntervalValue{Int64,Int64}\n(4,6) => 2\n IntervalTrees.IntervalValue{Int64,Int64}\n(5,7) => 2\n IntervalTrees.IntervalValue{Int64,Int64}\n(6,8) => 3"
},

{
    "location": "sequences/bioseq.html#The-pairwise_count-function-1",
    "page": "BioSequence",
    "title": "The pairwise_count function",
    "category": "section",
    "text": "Counting can also be done on a set of sequences in a pairwise manner with the count_pairwise function:julia> count_pairwise(Match, dna\"ATCGCCA-\", dna\"ATCGCCTA\", dna\"ATCGCCT-\", dna\"GTCGCCTA\")\n4×4 Array{Int64,2}:\n 0  6  7  5\n 6  0  7  7\n 7  7  0  6\n 5  7  6  0"
},

{
    "location": "sequences/bioseq.html#Iteration-1",
    "page": "BioSequence",
    "title": "Iteration",
    "category": "section",
    "text": "Sequences also work as iterators over symbols:julia> n = 0\n0\n\njulia> for nt in dna\"ATNGNNT\"\n           if nt == DNA_N\n               n += 1\n           end\n       end\n\njulia> n\n3\n"
},

{
    "location": "sequences/bioseq.html#Using-a-more-compact-sequence-representation-1",
    "page": "BioSequence",
    "title": "Using a more compact sequence representation",
    "category": "section",
    "text": "As we saw above, DNA and RNA sequences can store any ambiguous nucleotides like 'N'.  If you are sure that nucleotide sequences store unambiguous nucleotides only, you can save the memory space of sequences. DNAAlphabet{2} is an alphabet that uses two bits per base and limits to only unambiguous nucleotide symbols (ACGT in DNA and ACGU in RNA). To create a sequence of this alphabet, you need to explicitly pass DNAAlphabet{2} to BioSequence as its type parameter:julia> seq = BioSequence{DNAAlphabet{2}}(\"ACGT\")\n4nt DNA Sequence:\nACGT\nRecall that DNASequence is a type alias of BioSequence{DNAAlphabet{4}}, which uses four bits per base. That is, BioSequence{DNAAlphabet{2}} saves half of memory footprint compared to BioSequence{DNAAlphabet{4}}. If you need to handle reference genomes that are composed of five nucleotides, ACGTN, consider to use the ReferenceSequence type described in the Reference sequences section."
},

{
    "location": "sequences/bioseq.html#Defining-a-new-alphabet-1",
    "page": "BioSequence",
    "title": "Defining a new alphabet",
    "category": "section",
    "text": "The alphabet type parameter A of BioSequence{A} enables a user to extend functionality of BioSequence with minimum effort. As an example, definition of a new alphabet type representing a sequence of boolean values is shown below:julia> immutable BoolAlphabet <: Alphabet end\n\njulia> BioSequences.bitsof(::Type{BoolAlphabet}) = 1\n\njulia> BioSequences.eltype(::Type{BoolAlphabet}) = Bool\n\njulia> BioSequences.alphabet(::Type{BoolAlphabet}) = false:true\n\njulia> function BioSequences.encode(::Type{BoolAlphabet}, x::Bool)\n           return UInt64(ifelse(x, 0x01, 0x00))\n       end\n\njulia> function BioSequences.decode(::Type{BoolAlphabet}, x::UInt64)\n           if x > 0x01\n               throw(BioSequences.DecodeError(BoolAlphabet, x))\n           end\n           return ifelse(x == 0x00, false, true)\n       end\n"
},

{
    "location": "sequences/refseq.html#",
    "page": "Reference Sequences",
    "title": "Reference Sequences",
    "category": "page",
    "text": ""
},

{
    "location": "sequences/refseq.html#Reference-sequences-1",
    "page": "Reference Sequences",
    "title": "Reference sequences",
    "category": "section",
    "text": "DNASequence (alias of BioSequence{DNAAlphabet{4}}) is a flexible data structure but always consumes 4 bits per base, which will waste a large part of the memory space when storing reference genome sequences.  In such a case, ReferenceSequence is helpful because it compresses positions of 'N' symbols so that long DNA sequences are stored with almost 2 bits per base. An important limitation is that the ReferenceSequence type is immutable due to the compression. Other sequence-like operations are supported:julia> seq = ReferenceSequence(dna\"NNCGTATTTTCN\")\n12nt Reference Sequence:\nNNCGTATTTTCN\n\njulia> seq[1]\nDNA_N\n\njulia> seq[5]\nDNA_T\n\njulia> seq[2:6]\n5nt Reference Sequence:\nNCGTA\n\njulia> ReferenceSequence(dna\"ATGM\")  # DNA_M is not accepted\nERROR: ArgumentError: invalid symbol M ∉ {A,C,G,T,N} at 4\n in convert at /Users/kenta/.julia/v0.4/Bio/src/seq/refseq.jl:58\n in call at essentials.jl:56\n"
},

{
    "location": "sequences/kmer.html#",
    "page": "Nucleic acid k-mers",
    "title": "Nucleic acid k-mers",
    "category": "page",
    "text": "CurrentModule = BioSequences\nDocTestSetup = quote\n    using BioSequences\nend"
},

{
    "location": "sequences/kmer.html#BioSequences.each",
    "page": "Nucleic acid k-mers",
    "title": "BioSequences.each",
    "category": "Function",
    "text": "each(::Type{Kmer{T,k}}, seq::Sequence[, step=1])\n\nInitialize an iterator over all k-mers in a sequence seq skipping ambiguous nucleotides without changing the reading frame.\n\nArguments\n\nKmer{T,k}: k-mer type to enumerate.\nseq: a nucleotide sequence.\nstep=1: the number of positions between iterated k-mers\n\nExamples\n\n# iterate over DNA codons\nfor (pos, codon) in each(DNAKmer{3}, dna\"ATCCTANAGNTACT\", 3)\n    @show pos, codon\nend\n\n\n\n"
},

{
    "location": "sequences/kmer.html#BioSequences.canonical",
    "page": "Nucleic acid k-mers",
    "title": "BioSequences.canonical",
    "category": "Function",
    "text": "canonical(kmer::Kmer)\n\nReturn the canonical k-mer of x.\n\nA canonical k-mer is the numerical lesser of a k-mer and its reverse complement. This is useful in hashing/counting k-mers in data that is not strand specific, and thus observing k-mer is equivalent to observing its reverse complement.\n\n\n\n"
},

{
    "location": "sequences/kmer.html#BioSequences.neighbors",
    "page": "Nucleic acid k-mers",
    "title": "BioSequences.neighbors",
    "category": "Function",
    "text": "neighbors(kmer::Kmer)\n\nReturn an iterator through k-mers neighboring kmer on a de Bruijn graph.\n\n\n\n"
},

{
    "location": "sequences/kmer.html#Nucleic-acid-k-mers-1",
    "page": "Nucleic acid k-mers",
    "title": "Nucleic acid k-mers",
    "category": "section",
    "text": "A common strategy to simplify the analysis of sequence data is to operate or short k-mers, for size fixed size k. These can be packed into machine integers allowing extremely efficient code. The BioSequences module has built in support for representing short sequences in 64-bit integers. Besides being fixed length, Kmer types, unlike other sequence types cannot contain ambiguous symbols like 'N'.The Kmer{T,k} type parameterized on symbol type (T, either DNA, or RNA) and size k. For ease of writing code, two type aliases for each nucleotide type are defined and named as DNAKmer{k} and RNAKmer{k}:julia> DNAKmer(\"ACGT\")  # create a DNA 4-mer from a string\nDNA 4-mer:\nACGT\n\njulia> RNAKmer(\"ACGU\")  # create an RNA 4-mer from a string\nRNA 4-mer:\nACGU\n\njulia> kmer\"ACGT\" # DNA k-mers may also be written as literals\nDNA 4-mer:\nACGT\n\njulia> typeof(DNAKmer(\"ACGT\"))\nBioSequences.Kmer{BioSymbols.DNA,4}each\ncanonical\nneighbors"
},

{
    "location": "io/fasta.html#",
    "page": "FASTA formatted files",
    "title": "FASTA formatted files",
    "category": "page",
    "text": "CurrentModule = BioSequences\nDocTestSetup = quote\n    using BioSequences\nend"
},

{
    "location": "io/fasta.html#IO-FASTA-formatted-files-1",
    "page": "FASTA formatted files",
    "title": "IO - FASTA formatted files",
    "category": "section",
    "text": "FASTA is a text-based file format for representing biological sequences. A FASTA file stores a list of sequence records with name, description, and sequence.The template of a sequence record is:>{name} {description}?\n{sequence}Here is an example of a chromosomal sequence:>chrI chromosome 1\nCCACACCACACCCACACACCCACACACCACACCACACACCACACCACACC\nCACACACACACATCCTAACACTACCCTAACACAGCCCTAATCTAACCCTG"
},

{
    "location": "io/fasta.html#BioSequences.FASTA.Reader",
    "page": "FASTA formatted files",
    "title": "BioSequences.FASTA.Reader",
    "category": "Type",
    "text": "FASTA.Reader(input::IO; index=nothing)\n\nCreate a data reader of the FASTA file format.\n\nArguments\n\ninput: data source\nindex=nothing: filepath to a random access index (currently fai is supported)\n\n\n\n"
},

{
    "location": "io/fasta.html#BioSequences.FASTA.Writer",
    "page": "FASTA formatted files",
    "title": "BioSequences.FASTA.Writer",
    "category": "Type",
    "text": "FASTA.Writer(output::IO; width=70)\n\nCreate a data writer of the FASTA file format.\n\nArguments\n\noutput: data sink\nwidth=70: wrapping width of sequence characters\n\n\n\n"
},

{
    "location": "io/fasta.html#BioSequences.FASTA.Record",
    "page": "FASTA formatted files",
    "title": "BioSequences.FASTA.Record",
    "category": "Type",
    "text": "FASTA.Record()\n\nCreate an unfilled FASTA record.\n\n\n\nFASTA.Record(data::Vector{UInt8})\n\nCreate a FASTA record object from data.\n\nThis function verifies and indexes fields for accessors. Note that the ownership of data is transferred to a new record object.\n\n\n\nFASTA.Record(str::AbstractString)\n\nCreate a FASTA record object from str.\n\nThis function verifies and indexes fields for accessors.\n\n\n\nFASTA.Record(identifier, sequence)\n\nCreate a FASTA record object from identifier and sequence.\n\n\n\nFASTA.Record(identifier, description, sequence)\n\nCreate a FASTA record object from identifier, description and sequence.\n\n\n\n"
},

{
    "location": "io/fasta.html#BioSequences.FASTA.hasidentifier",
    "page": "FASTA formatted files",
    "title": "BioSequences.FASTA.hasidentifier",
    "category": "Function",
    "text": "hasidentifier(record::Record)\n\nChecks whether or not the record has an identifier.\n\n\n\n"
},

{
    "location": "io/fasta.html#BioSequences.FASTA.identifier",
    "page": "FASTA formatted files",
    "title": "BioSequences.FASTA.identifier",
    "category": "Function",
    "text": "identifier(record::Record)::String\n\nGet the sequence identifier of record.\n\n\n\n"
},

{
    "location": "io/fasta.html#BioSequences.FASTA.hasdescription",
    "page": "FASTA formatted files",
    "title": "BioSequences.FASTA.hasdescription",
    "category": "Function",
    "text": "hasdescription(record::Record)\n\nChecks whether or not the record has a description.\n\n\n\n"
},

{
    "location": "io/fasta.html#BioSequences.FASTA.description",
    "page": "FASTA formatted files",
    "title": "BioSequences.FASTA.description",
    "category": "Function",
    "text": "description(record::Record)::String\n\nGet the description of record.\n\n\n\n"
},

{
    "location": "io/fasta.html#BioSequences.FASTA.hassequence",
    "page": "FASTA formatted files",
    "title": "BioSequences.FASTA.hassequence",
    "category": "Function",
    "text": "hassequence(record::Record)\n\nChecks whether or not a sequence record contains a sequence.\n\n\n\n"
},

{
    "location": "io/fasta.html#BioSequences.FASTA.sequence-Tuple{BioSequences.FASTA.Record,UnitRange{Int64}}",
    "page": "FASTA formatted files",
    "title": "BioSequences.FASTA.sequence",
    "category": "Method",
    "text": "sequence(record::Record, [part::UnitRange{Int}])\n\nGet the sequence of record.\n\nThis function infers the sequence type from the data. When it is wrong or unreliable, use sequence(::Type{S}, record::Record).  If part argument is given, it returns the specified part of the sequence.\n\n\n\n"
},

{
    "location": "io/fasta.html#Readers-and-Writers-1",
    "page": "FASTA formatted files",
    "title": "Readers and Writers",
    "category": "section",
    "text": "The reader and writer for FASTA formatted files, are found within the BioSequences.FASTA module.FASTA.Reader\nFASTA.WriterThey can be created with IOStreams:r = FASTA.Reader(open(\"MyInput.fasta\", \"r\"))\nw = FASTA.Writer(open(\"MyFile.fasta\", \"w\"))Usually sequence records will be read sequentially from a file by iteration.using BioSequences\nreader = FASTA.Reader(open(\"hg38.fa\", \"r\"))\nfor record in reader\n    # Do something\nend\nclose(reader)But if the FASTA file has an auxiliary index file formatted in fai, the reader supports random access to FASTA records, which would be useful when accessing specific parts of a huge genome sequence:reader = open(FASTAReader, \"sacCer.fa\", index=\"sacCer.fa.fai\")\nchrIV = reader[\"chrIV\"]  # directly read sequences called chrIV.Reading in a sequence from a FASTA formatted file will give you a variable of type FASTA.Record.FASTA.RecordVarious getters and setters are available for FASTA.Records:FASTA.hasidentifier\nFASTA.identifier\nFASTA.hasdescription\nFASTA.description\nFASTA.hassequence\nFASTA.sequence(record::FASTA.Record, [part::UnitRange{Int}])To write a BioSequence to FASTA file, you first have to create a FASTA.Record:using BioSequences\nx = dna\"aaaaatttttcccccggggg\"\nrec = FASTA.Record(\"MySeq\", x)\nw = FASTA.Writer(open(\"MyFile.fasta\", \"w\"))\nwrite(w, rec)As always with julia IO types, remember to close your file readers and writer after you are finished."
},

{
    "location": "io/fastq.html#",
    "page": "FASTQ formatted files",
    "title": "FASTQ formatted files",
    "category": "page",
    "text": "CurrentModule = BioSequences\nDocTestSetup = quote\n    using BioSequences\nend"
},

{
    "location": "io/fastq.html#IO-FASTQ-formatted-files-1",
    "page": "FASTQ formatted files",
    "title": "IO - FASTQ formatted files",
    "category": "section",
    "text": "FASTQ is a text-based file format for representing DNA sequences along with qualities for each base. A FASTQ file stores a list of sequence records in the following format:@{name} {description}?\n{sequence}\n+\n{qualities}Here is an example of one record from a FASTQ file:@FSRRS4401BE7HA\ntcagTTAAGATGGGAT\n+\n###EEEEEEEEE##E#"
},

{
    "location": "io/fastq.html#BioSequences.FASTQ.Reader",
    "page": "FASTQ formatted files",
    "title": "BioSequences.FASTQ.Reader",
    "category": "Type",
    "text": "FASTQ.Reader(input::IO; fill_ambiguous=nothing)\n\nCreate a data reader of the FASTQ file format.\n\nArguments\n\ninput: data source\nfill_ambiguous=nothing: fill ambiguous symbols with the given symbol\n\n\n\n"
},

{
    "location": "io/fastq.html#BioSequences.FASTQ.Writer",
    "page": "FASTQ formatted files",
    "title": "BioSequences.FASTQ.Writer",
    "category": "Type",
    "text": "FASTQ.Writer(output::IO; quality_header=false)\n\nCreate a data writer of the FASTQ file format.\n\nArguments\n\noutput: data sink\nquality_header=false: output the title line at the third line just after '+'\n\n\n\n"
},

{
    "location": "io/fastq.html#BioSequences.FASTQ.Record",
    "page": "FASTQ formatted files",
    "title": "BioSequences.FASTQ.Record",
    "category": "Type",
    "text": "FASTQ.Record()\n\nCreate an unfilled FASTQ record.\n\n\n\nFASTQ.Record(data::Vector{UInt8})\n\nCreate a FASTQ record object from data.\n\nThis function verifies and indexes fields for accessors. Note that the ownership of data is transferred to a new record object.\n\n\n\nFASTQ.Record(str::AbstractString)\n\nCreate a FASTQ record object from str.\n\nThis function verifies and indexes fields for accessors.\n\n\n\nFASTQ.Record(identifier, sequence, quality; offset=33)\n\nCreate a FASTQ record from identifier, sequence and quality.\n\n\n\nFASTQ.Record(identifier, description, sequence, quality; offset=33)\n\nCreate a FASTQ record from identifier, description, sequence and quality.\n\n\n\n"
},

{
    "location": "io/fastq.html#BioSequences.FASTQ.hasidentifier",
    "page": "FASTQ formatted files",
    "title": "BioSequences.FASTQ.hasidentifier",
    "category": "Function",
    "text": "hasidentifier(record::Record)\n\nChecks whether or not the record has an identifier.\n\n\n\n"
},

{
    "location": "io/fastq.html#BioSequences.FASTQ.identifier",
    "page": "FASTQ formatted files",
    "title": "BioSequences.FASTQ.identifier",
    "category": "Function",
    "text": "identifier(record::Record)::String\n\nGet the sequence identifier of record.\n\n\n\n"
},

{
    "location": "io/fastq.html#BioSequences.FASTQ.hasdescription",
    "page": "FASTQ formatted files",
    "title": "BioSequences.FASTQ.hasdescription",
    "category": "Function",
    "text": "hasdescription(record::Record)\n\nChecks whether or not the record has a description.\n\n\n\n"
},

{
    "location": "io/fastq.html#BioSequences.FASTQ.description",
    "page": "FASTQ formatted files",
    "title": "BioSequences.FASTQ.description",
    "category": "Function",
    "text": "description(record::Record)::String\n\nGet the description of record.\n\n\n\n"
},

{
    "location": "io/fastq.html#BioSequences.FASTQ.hassequence",
    "page": "FASTQ formatted files",
    "title": "BioSequences.FASTQ.hassequence",
    "category": "Function",
    "text": "hassequence(record::Record)\n\nChecks whether or not a sequence record contains a sequence.\n\n\n\n"
},

{
    "location": "io/fastq.html#BioSequences.FASTQ.sequence-Tuple{BioSequences.FASTQ.Record,UnitRange{Int64}}",
    "page": "FASTQ formatted files",
    "title": "BioSequences.FASTQ.sequence",
    "category": "Method",
    "text": "sequence(record::Record, [part::UnitRange{Int}])::BioSequences.DNASequence\n\nGet the sequence of record.\n\n\n\n"
},

{
    "location": "io/fastq.html#BioSequences.FASTQ.hasquality",
    "page": "FASTQ formatted files",
    "title": "BioSequences.FASTQ.hasquality",
    "category": "Function",
    "text": "hasquality(record::Record)\n\nCheck whether the given FASTQ record has a quality string.\n\n\n\n"
},

{
    "location": "io/fastq.html#BioSequences.FASTQ.quality",
    "page": "FASTQ formatted files",
    "title": "BioSequences.FASTQ.quality",
    "category": "Function",
    "text": "quality(record::Record, [offset::Integer=33, [part::UnitRange]])::Vector{UInt8}\n\nGet the base quality of record.\n\n\n\nquality(record::Record, encoding_name::Symbol, [part::UnitRange])::Vector{UInt8}\n\nGet the base quality of record by decoding with encoding_name.\n\nThe encoding_name can be either :sanger, :solexa, :illumina13, :illumina15, or :illumina18.\n\n\n\n"
},

{
    "location": "io/fastq.html#BioSequences.FASTQ.Record-Tuple{AbstractString,Union{AbstractString, Void},Any,Array{T,1} where T}",
    "page": "FASTQ formatted files",
    "title": "BioSequences.FASTQ.Record",
    "category": "Method",
    "text": "FASTQ.Record(identifier, description, sequence, quality; offset=33)\n\nCreate a FASTQ record from identifier, description, sequence and quality.\n\n\n\n"
},

{
    "location": "io/fastq.html#Readers-and-Writers-1",
    "page": "FASTQ formatted files",
    "title": "Readers and Writers",
    "category": "section",
    "text": "The reader and writer for FASTQ formatted files, are found within the BioSequences.FASTQ module.FASTQ.Reader\nFASTQ.WriterThey can be created with IOStreams:r = FASTQ.Reader(open(\"MyInput.fastq\", \"r\"))\nw = FASTQ.Writer(open(\"MyFile.fastq\", \"w\"))Note that FASTQ.Reader does not support line-wraps within sequence and quality. Usually sequence records will be read sequentially from a file by iteration.using BioSequences\nreader = FASTQ.Reader(open(\"hg38.fastq\", \"r\"))\nfor record in reader\n    # Do something\nend\nclose(reader)Reading in a record from a FASTQ formatted file will give you a variable of type FASTQ.Record.FASTQ.RecordVarious getters and setters are available for FASTQ.Records:FASTQ.hasidentifier\nFASTQ.identifier\nFASTQ.hasdescription\nFASTQ.description\nFASTQ.hassequence\nFASTQ.sequence(record::FASTQ.Record, [part::UnitRange{Int}])\nFASTQ.hasquality\nFASTQ.qualityTo write a BioSequence to FASTQ file, you first have to create a FASTQ.Record:FASTQ.Record(identifier::AbstractString, description::Union{AbstractString,Void}, sequence, quality::Vector; offset=33)As always with julia IO types, remember to close your file readers and writer after you are finished."
},

{
    "location": "io/twobit.html#",
    "page": "2bit formatted files",
    "title": "2bit formatted files",
    "category": "page",
    "text": "CurrentModule = BioSequences\nDocTestSetup = quote\n    using BioSequences\nend"
},

{
    "location": "io/twobit.html#IO-2bit-formatted-files-1",
    "page": "2bit formatted files",
    "title": "IO - 2bit formatted files",
    "category": "section",
    "text": "2bit is a binary file format designed for storing a genome consists of multiple chromosomal sequences. The reading speed is often an order of magnitude faster than that of FASTA and the file size is smaller. However, since the .2bit file format is specialized for genomic sequences, it cannot store either RNA or amino acid sequences."
},

{
    "location": "io/twobit.html#BioSequences.TwoBit.Reader",
    "page": "2bit formatted files",
    "title": "BioSequences.TwoBit.Reader",
    "category": "Type",
    "text": "TwoBit.Reader(input::IO)\n\nCreate a data reader of the 2bit file format.\n\nArguments\n\ninput: data source\n\n\n\n"
},

{
    "location": "io/twobit.html#BioSequences.TwoBit.Writer",
    "page": "2bit formatted files",
    "title": "BioSequences.TwoBit.Writer",
    "category": "Type",
    "text": "TwoBitWriter(output::IO, names::AbstractVector)\n\nCreate a data writer of the 2bit file format.\n\nArguments\n\noutput: data sink\nnames: a vector of sequence names written to output\n\n\n\n"
},

{
    "location": "io/twobit.html#BioSequences.TwoBit.Record",
    "page": "2bit formatted files",
    "title": "BioSequences.TwoBit.Record",
    "category": "Type",
    "text": "TwoBit.Record()\n\nCreate an unfilled 2bit record.\n\n\n\nRecord()\n\nPrepare a record for writing to a 2bit formatted file.\n\nNeeds a name, a sequence, and (optionally) masks: a vector of ranges that delineate masked regions of sequence.\n\n\n\n"
},

{
    "location": "io/twobit.html#BioSequences.TwoBit.Record",
    "page": "2bit formatted files",
    "title": "BioSequences.TwoBit.Record",
    "category": "Type",
    "text": "Record()\n\nPrepare a record for writing to a 2bit formatted file.\n\nNeeds a name, a sequence, and (optionally) masks: a vector of ranges that delineate masked regions of sequence.\n\n\n\n"
},

{
    "location": "io/twobit.html#Readers-and-Writers-1",
    "page": "2bit formatted files",
    "title": "Readers and Writers",
    "category": "section",
    "text": "The reader and writer for 2bit formatted files, are found within the BioSequences.TwoBit module.TwoBit.Reader\nTwoBit.WriterThe 2bit reader supports random access using an index included in the header section of a .2bit file:reader = TwoBit.Reader(open(\"sacCer.2bit\", \"r\"))\nchrIV = reader[\"chrIV\"] # directly read chromosome 4If you want to know the names of the sequences available in the file, you can use the seqnames method on the reader.seqnames(reader)Reading from a TwoBit.Reader will yield a TwoBit.Record type variable:TwoBit.RecordTo write a sequence to a TwoBit file, first a record must be created.TwoBit.Record(name::AbstractString, seq::BioSequences.Sequence, masks = Nullable{Vector{UnitRange{Int}}}())"
},

{
    "location": "search.html#",
    "page": "Searching",
    "title": "Searching",
    "category": "page",
    "text": "CurrentModule = BioSequences\nDocTestSetup = quote\n    using BioSequences\nend"
},

{
    "location": "search.html#Sequence-search-1",
    "page": "Searching",
    "title": "Sequence search",
    "category": "section",
    "text": "Three kinds of on-line search functions are provided:Exact search\nApproximate search\nRegular expression searchThese are all specialized for biological sequences and ambiguities of symbols are considered."
},

{
    "location": "search.html#Exact-search-1",
    "page": "Searching",
    "title": "Exact search",
    "category": "section",
    "text": "Exact search functions search for an occurrence of the query symbol or sequence. Four functions, search, searchindex, rsearch, and rsearchindex are available:julia> seq = dna\"ACAGCGTAGCT\";\n\njulia> search(seq, DNA_G)  # search a query symbol\n4:4\n\njulia> query = dna\"AGC\";\n\njulia> search(seq, query)  # search a query sequence\n3:5\n\njulia> searchindex(seq, query)\n3\n\njulia> rsearch(seq, query)  # similar to `search` but in the reverse direction\n8:10\n\njulia> rsearchindex(seq, query)  # similar to `searchindex` but in the reverse direction\n8\nThese search functions take ambiguous symbols into account. That is, if two symbols are compatible (e.g. DNA_A and DNA_N), they match when searching an occurrence. In the following example, 'N' is a wild card that matches any symbols:julia> search(dna\"ACNT\", DNA_N)  # 'A' matches 'N'\n1:1\n\njulia> search(dna\"ACNT\", dna\"CGT\")  # 'N' matches 'G'\n2:4\n\njulia> search(dna\"ACGT\", dna\"CNT\")  # 'G' matches 'N'\n2:4\nThe exact sequence search needs preprocessing phase of query sequence before searching phase. This would be enough fast for most search applications. But when searching a query sequence to large amounts of target sequences, caching the result of preprocessing may save time. The ExactSearchQuery creates such a preprocessed query object and is applicable to the search functions:julia> query = ExactSearchQuery(dna\"ATT\");\n\njulia> search(dna\"ATTTATT\", query)\n1:3\n\njulia> rsearch(dna\"ATTTATT\", query)\n5:7\n"
},

{
    "location": "search.html#Approximate-search-1",
    "page": "Searching",
    "title": "Approximate search",
    "category": "section",
    "text": "The approximate search is similar to the exact search but allows a specific number of errors. That is, it tries to find a subsequence of the target sequence within a specific Levenshtein distance of the query sequence:julia> seq = dna\"ACAGCGTAGCT\";\n\njulia> approxsearch(seq, dna\"AGGG\", 0)  # nothing matches with no errors\n0:-1\n\njulia> approxsearch(seq, dna\"AGGG\", 1)  # seq[3:5] matches with one error\n3:6\n\njulia> approxsearch(seq, dna\"AGGG\", 2)  # seq[1:4] matches with two errors\n1:4\nLike the exact search functions, four kinds of functions (approxsearch, approxsearchindex, approxrsearch, and approxrsearchindex) are available:julia> seq = dna\"ACAGCGTAGCT\"; pat = dna\"AGGG\";\n\njulia> approxsearch(seq, pat, 2)        # return the range (forward)\n1:4\n\njulia> approxsearchindex(seq, pat, 2)   # return the starting index (forward)\n1\n\njulia> approxrsearch(seq, pat, 2)       # return the range (backward)\n8:11\n\njulia> approxrsearchindex(seq, pat, 2)  # return the starting index (backward)\n8\nPreprocessing can be cached in an ApproximateSearchQuery object:julia> query = ApproximateSearchQuery(dna\"AGGG\");\n\njulia> approxsearch(dna\"AAGAGG\", query, 1)\n2:5\n\njulia> approxsearch(dna\"ACTACGT\", query, 2)\n4:6\n"
},

{
    "location": "search.html#Regular-expression-search-1",
    "page": "Searching",
    "title": "Regular expression search",
    "category": "section",
    "text": "Query patterns can be described in regular expressions. The syntax supports a subset of Perl and PROSITE's notation.The Perl-like syntax starts with biore (biological regular expression) and ends with a symbol option: \"dna\", \"rna\" or \"aa\". For example, biore\"A+\"dna is a regular expression for DNA sequences and biore\"A+\"aa is for amino acid sequences. The symbol options can be abbreviated to its first character: \"d\", \"r\" or \"a\", respectively.Here are examples of using the regular expression for BioSequences:julia> match(biore\"A+C*\"dna, dna\"AAAACC\")\nNullable{BioSequences.RE.RegexMatch{BioSequences.BioSequence{BioSequences.DNAAlphabet{4}}}}(RegexMatch(\"AAAACC\"))\n\njulia> match(biore\"A+C*\"d, dna\"AAAACC\")\nNullable{BioSequences.RE.RegexMatch{BioSequences.BioSequence{BioSequences.DNAAlphabet{4}}}}(RegexMatch(\"AAAACC\"))\n\njulia> ismatch(biore\"A+C*\"dna, dna\"AAC\")\ntrue\n\njulia> ismatch(biore\"A+C*\"dna, dna\"C\")\nfalse\nmatch always returns a Nullable object and it should be null if no match is found.The table below summarizes available syntax elements.Syntax Description Example\n| alternation \"A|T\" matches \"A\" and \"T\"\n* zero or more times repeat \"TA*\" matches \"T\", \"TA\" and \"TAA\"\n+ one or more times repeat \"TA+\" matches \"TA\" and \"TAA\"\n? zero or one time \"TA?\" matches \"T\" and \"TA\"\n{n,} n or more times repeat \"A{3,}\" matches \"AAA\" and \"AAAA\"\n{n,m} n-m times repeat \"A{3,5}\" matches \"AAA\", \"AAAA\" and \"AAAAA\"\n^ the start of the sequence \"^TAN*\" matches \"TATGT\"\n$ the end of the sequence \"N*TA$\" matches \"GCTA\"\n(...) pattern grouping \"(TA)+\" matches \"TA\" and \"TATA\"\n[...] one of symbols \"[ACG]+\" matches \"AGGC\"eachmatch, matchall, and search are also defined like usual strings:julia> matchall(biore\"TATA*?\"d, dna\"TATTATAATTA\")  # overlap (default)\n4-element Array{BioSequences.BioSequence{BioSequences.DNAAlphabet{4}},1}:\n TAT  \n TAT  \n TATA\n TATAA\n\njulia> matchall(biore\"TATA*\"d, dna\"TATTATAATTA\", false)  # no overlap\n2-element Array{BioSequences.BioSequence{BioSequences.DNAAlphabet{4}},1}:\n TAT  \n TATAA\n\njulia> search(dna\"TATTATAATTA\", biore\"TATA*\"d)\n1:3\n\njulia> search(dna\"TATTATAATTA\", biore\"TATA*\"d, 2)\n4:8\nNotewothy differences from strings are:Ambiguous characters match any compatible characters (e.g. biore\"N\"d is equivalent to biore\"[ACGT]\"d).\nWhitespaces are ignored (e.g. biore\"A C G\"d is equivalent to biore\"ACG\"d).The PROSITE notation is described in ScanProsite - user manual. The syntax supports almost all notations including the extended syntax. The PROSITE notation starts with prosite prefix and no symbol option is needed because it always describes patterns of amino acid sequences:julia> match(prosite\"[AC]-x-V-x(4)-{ED}\", aa\"CPVPQARG\")\nNullable{BioSequences.RE.RegexMatch{BioSequences.BioSequence{BioSequences.AminoAcidAlphabet}}}(RegexMatch(\"CPVPQARG\"))\n\njulia> match(prosite\"[AC]xVx(4){ED}\", aa\"CPVPQARG\")\nNullable{BioSequences.RE.RegexMatch{BioSequences.BioSequence{BioSequences.AminoAcidAlphabet}}}(RegexMatch(\"CPVPQARG\"))\n"
},

{
    "location": "search.html#Position-weight-matrix-search-1",
    "page": "Searching",
    "title": "Position weight matrix search",
    "category": "section",
    "text": "A motif can also be specified using position weight matrix (PWM) in a probabilistic way. search(seq, pwm, threshold) method searches for the first position in the sequence where a score calculated using the PWM is greater than or equal to the threshold. More formally, denoting the sequence as S and the PWM value of symbol s at position j as M_sj, the score starting from a position p is defined asoperatornamescore(S p) = sum_i=1^L M_Sp+i-1iand search(S, M, t) returns the smallest p that satisfies operatornamescore(S p) ge t.There are two kinds of matrices in this package: PFM and PWM. The PFM type is a position frequency matrix and stores symbol frequencies for each position. The PWM is a position weight matrix and stores symbol scores for each position. You can create a PFM from a set of sequences with the same length and then create a PWM from the PFM object.julia> kmers = DNAKmer.([\"TTA\", \"CTA\", \"ACA\", \"TCA\", \"GTA\"])\n5-element Array{BioSequences.Kmer{BioSymbols.DNA,3},1}:\n TTA\n CTA\n ACA\n TCA\n GTA\n\njulia> pfm = PFM(kmers)  # sequence set => PFM\n4×3 BioSequences.PFM{BioSymbols.DNA,Int64}:\n A  1  0  5\n C  1  2  0\n G  1  0  0\n T  2  3  0\n\njulia> pwm = PWM(pfm)  # PFM => PWM\n4×3 BioSequences.PWM{BioSymbols.DNA,Float64}:\n A -0.321928 -Inf       2.0\n C -0.321928  0.678072 -Inf\n G -0.321928 -Inf      -Inf\n T  0.678072  1.26303  -Inf\n\njulia> pwm = PWM(pfm .+ 0.01)  # add pseudo counts to avoid infinite values\n4×3 BioSequences.PWM{BioSymbols.DNA,Float64}:\n A -0.319068 -6.97728   1.99139\n C -0.319068  0.673772 -6.97728\n G -0.319068 -6.97728  -6.97728\n T  0.673772  1.25634  -6.97728\n\njulia> pwm = PWM(pfm .+ 0.01, prior=[0.2, 0.3, 0.3, 0.2])  # GC-rich prior\n4×3 BioSequences.PWM{BioSymbols.DNA,Float64}:\n A  0.00285965 -6.65535   2.31331\n C -0.582103    0.410737 -7.24031\n G -0.582103   -7.24031  -7.24031\n T  0.9957      1.57827  -6.65535\nThe PWM_sj matrix is computed from PFM_sj and the prior probability p(s) as follows ([Wasserman2004]):beginalign\n    PWM_sj = log_2 fracp(sj)p(s) \n    p(sj)  = fracPFM_sjsum_s PFM_sj\nendalign[Wasserman2004]: https://doi.org/10.1038/nrg1315"
},

{
    "location": "composition.html#",
    "page": "Sequence Composition",
    "title": "Sequence Composition",
    "category": "page",
    "text": "CurrentModule = BioSequences\nDocTestSetup = quote\n    using BioSequences\nend"
},

{
    "location": "composition.html#BioSequences.Composition",
    "page": "Sequence Composition",
    "title": "BioSequences.Composition",
    "category": "Type",
    "text": "Sequence composition.\n\nThis is a subtype of Associative{T,Int}, and the getindex method returns the number of occurrences of a symbol or a k-mer.\n\n\n\n"
},

{
    "location": "composition.html#BioSequences.composition",
    "page": "Sequence Composition",
    "title": "BioSequences.composition",
    "category": "Function",
    "text": "composition(seq | kmer_iter)\n\nCalculate composition of biological symbols in seq or k-mers in kmer_iter.\n\n\n\ncomposition(iter)\n\nA generalised composition algorithm, which computes the number of unique items produced by an iterable.\n\nExample\n\n\n# Example, counting unique sequences.\n\njulia> a = dna\"AAAAAAAATTTTTT\"\n14nt DNA Sequence:\nAAAAAAAATTTTTT\n\njulia> b = dna\"AAAAAAAATTTTTT\"\n14nt DNA Sequence:\nAAAAAAAATTTTTT\n\njulia> c = a[5:10]\n6nt DNA Sequence:\nAAAATT\n\njulia> composition([a, b, c])\nVector{BioSequences.BioSequence{BioSequences.DNAAlphabet{4}}} Composition:\n  AAAATT         => 1\n  AAAAAAAATTTTTT => 2\n\n\n\n"
},

{
    "location": "composition.html#Sequence-composition-1",
    "page": "Sequence Composition",
    "title": "Sequence composition",
    "category": "section",
    "text": "There are many instances in analyzing sequence data where you will want to know about the composition of your sequences.For example, for a given sequence, you may want to count how many of each possible Kmer, is present in the sequence. This would be important if - for instance - you wanted to analyze the Kmer spectra of your data. Alternatively you might have a collection of sequences, and may want to count how many of each unique sequence you have in your collection. This would be important if - for instance - your collection of sequences were from a population sample, and you wanted to compute the allele or genotype frequencies for the population.Whatever the application, BioSequences provides a method called composition, and a parametric struct called Composition to both compute, and handle the results of such sequence composition calculations.Composition{T}\ncompositionFor example to get the nucleotide composition of a sequence:julia> comp = composition(dna\"ACGAG\")\nDNA Composition:\n  DNA_A => 2\n  DNA_G => 2\n  DNA_C => 1\n\njulia> comp[DNA_A]\n2\n\njulia> comp[DNA_T]\n0\nComposition structs behave like an associative collection, such as a Dict. But there are a few differences:The getindex method for Composition structs is overloaded to return a default value of 0, if a key is used that is not present in the Composition.\nThe merge! method for two Composition structs adds counts together, unlike the merge! method for other associative containers, which would overwrite the counts.merge! is used to accumulate composition statistics of multiple sequences:# initiaize an empty composition counter\ncomp = composition(dna\"\");\n\n# iterate over sequences and accumulate composition statistics into `comp`\nfor seq in seqs\n    merge!(comp, composition(seq))\nend\n\n# or functional programming style in one line\nfoldl((x, y) -> merge(x, composition(y)), composition(dna\"\"), seqs)composition is also applicable to a k-mer iterator:julia> comp = composition(each(DNAKmer{4}, dna\"ACGT\"^100));\n\njulia> comp[DNAKmer(\"ACGT\")]\n100\n\njulia> comp[DNAKmer(\"CGTA\")]\n99\n"
},

{
    "location": "demultiplexer.html#",
    "page": "Demultiplexing",
    "title": "Demultiplexing",
    "category": "page",
    "text": "CurrentModule = BioSequences\nDocTestSetup = quote\n    using BioSequences\nend"
},

{
    "location": "demultiplexer.html#Sequence-demultiplexing-1",
    "page": "Demultiplexing",
    "title": "Sequence demultiplexing",
    "category": "section",
    "text": "Multiplex sequencing is a technology to sequence multiple samples at the same time on a high-throughput DNA sequencer. Samples are distinguished by the short prefix of a DNA sequence called DNA barcode. The BioSequences offers the Demultiplexer type and the demultiplex function to identify the DNA barcode of a longer DNA sequence allowing small errors.In the following example, four kinds of DNA sequences of length 4 are used as DNA barcodes. Demultiplexer takes these barcodes as its first argument with a few options:julia> barcodes = DNASequence[\"ATGG\", \"CAGA\", \"GGAA\", \"TACG\"];\n\njulia> dplxr = Demultiplexer(barcodes, n_max_errors=1, distance=:hamming)\nBioSequences.Demultiplexer{BioSequences.BioSequence{BioSequences.DNAAlphabet{4}}}:\n  distance: hamming\n  number of barcodes: 4\n  number of correctable errors: 1DocTestSetup = quote\n    using BioSequences\n    barcodes = DNASequence[\"ATGG\", \"CAGA\", \"GGAA\", \"TACG\"];\n    dplxr = Demultiplexer(barcodes, n_max_errors=1, distance=:hamming);\nendn_max_errors specifies the number of maximum correctable errors in a barcode. The type of correctable errors depends on the distance parameter. When distance = :hamming as shown above only substitutions are correctable. When distance = :levenshtein substitutions, deletions, and insertions are correctable. The user is responsible for keeping enough distances among barcodes; Demultiplexer will throw an exception if two barcodes are within n_max_errors * 2.The demultiplex function takes a demultiplexer object and a DNA sequence, and returns a tuple of a barcode index and a distance between the original barcode sequence and the prefix sequence:julia> demultiplex(dplxr, dna\"ATGGCGNT\")  # 1st barcode with no errors\n(1, 0)\n\njulia> demultiplex(dplxr, dna\"CAGGCGNT\")  # 2nd barcode with one error\n(2, 1)\n\njulia> demultiplex(dplxr, dna\"GGAACGNT\")  # 3rd barcode with no errors\n(3, 0)\n\njulia> demultiplex(dplxr, dna\"TGACCGNT\")  # no matching barcode\n(0, -1)\nThe optional third argument controls the search strategy. demultiplex uses an index to search the closest barcode within n_max_errors in the barcode set and returns it if any by default. If the third argument is true it falls back to a linear search after the index search and returns one of the closest barcodes at random. The next example shows the difference of these two strategies:julia> demultiplex(dplxr, dna\"TGACCGNT\", false)  # linear search off (default)\n(0, -1)\n\njulia> demultiplex(dplxr, dna\"TGACCGNT\", true)   # linear search on\n(3, 2)\n"
},

{
    "location": "contributing.html#",
    "page": "Contributing",
    "title": "Contributing",
    "category": "page",
    "text": ""
},

{
    "location": "contributing.html#Contributing-1",
    "page": "Contributing",
    "title": "Contributing",
    "category": "section",
    "text": "We appreciate contributions from users including reporting bugs, fixing issues, improving performance and adding new features.If you have a question about contributing or using this package, you are encouraged to use the Bio category of the Julia discourse site.Detailed guidance for contributing to all BioJulia packages is provided at the BioJulia Contribution Documentation.Here we list specific details about contributing and maintainership pertaining specifically to the BioSequences.jl package."
},

{
    "location": "contributing.html#Named-maintainers-1",
    "page": "Contributing",
    "title": "Named maintainers",
    "category": "section",
    "text": "The named maintainers of this package are Kenta Sato and Ben Ward. It is their responsibility to make final choices about pull requests and issues, although because of our community structure, you will find other maintainers assisting them."
},

{
    "location": "contributing.html#Branching-model-1",
    "page": "Contributing",
    "title": "Branching model",
    "category": "section",
    "text": "The branching model used to develop and make releases of this package is the OneFlow model summarized in the BioJulia Contribution Documentation"
},

]}
