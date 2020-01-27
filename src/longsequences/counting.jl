###
### LongSequence specific specializations of src/biosequence/counting.jl
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

# Counting GC positions
let
    @info "Compiling bit-parallel GC counter for LongSequence{<:NucleicAcidAlphabet}"
    counter = :(n += gc_bitcount(chunk, BitsPerSymbol(seq)))
    compile_bitpar(
        :count_gc_bitpar,
        arguments   = (:(seq::LongSequence{<:NucleicAcidAlphabet}),),
        init_code   = :(n = 0),
        head_code   = counter,
        body_code   = counter,
        tail_code   = counter,
        return_code = :(return n)
    ) |> eval
end

Base.count(::typeof(isGC), seq::LongSequence{<:NucleicAcidAlphabet}) = count_gc_bitpar(seq)

# Counting mismatches
let
    @info "Compiling bit-parallel mismatch counter for LongSequence{<:NucleicAcidAlphabet}"
    
    counter = :(count += mismatch_bitcount(x, y, A()))
    
    compile_2seq_bitpar(
        :count_mismatches_bitpar,
        arguments = (:(seqa::LongSequence{A}), :(seqb::LongSequence{A})),
        parameters = (:(A<:NucleicAcidAlphabet),),
        init_code = :(count = 0),
        head_code = counter,
        body_code = counter,
        tail_code = counter,
        return_code = :(return count)
    ) |> eval
end
Base.count(::typeof(!=), seqa::LongSequence{A}, seqb::LongSequence{A}) where {A<:NucleicAcidAlphabet} = count_mismatches_bitpar(seqa, seqb)
Base.count(::typeof(!=), seqa::LongNucleotideSequence, seqb::LongNucleotideSequence) = count(!=, promote(seqa, seqb)...)

# Counting matches
let
    @info "Compiling bit-parallel match counter for LongSequence{<:NucleicAcidAlphabet}"
    
    counter = :(count += match_bitcount(x, y, A()))
    
    count_empty = quote
        count += match_bitcount(x, y, A())
        nempty = div(64, bits_per_symbol(A())) - div(offs, bits_per_symbol(A()))
        count -= nempty
    end
    
    compile_2seq_bitpar(
        :count_matches_bitpar,
        arguments = (:(seqa::LongSequence{A}), :(seqb::LongSequence{A})),
        parameters = (:(A<:NucleicAcidAlphabet),),
        init_code = :(count = 0),
        head_code = count_empty,
        body_code = counter,
        tail_code = count_empty,
        return_code = :(return count)
    ) |> eval
end
Base.count(::typeof(==), seqa::LongSequence{A}, seqb::LongSequence{A}) where {A<:NucleicAcidAlphabet} = count_matches_bitpar(seqa, seqb)
Base.count(::typeof(==), seqa::LongNucleotideSequence, seqb::LongNucleotideSequence) = count(==, promote(seqa, seqb)...)

# Counting ambiguous sites
# ------------------------
let
    @info "Compiling bit-parallel ambiguity counter..."
    @info "\tFor a single LongSequence{<:NucleicAcidAlphabet}"
    
    counter = :(count += ambiguous_bitcount(chunk, Alphabet(seq)))
    
    compile_bitpar(
        :count_ambiguous_bitpar,
        arguments   = (:(seq::LongSequence{<:NucleicAcidAlphabet}),),
        init_code   = :(count = 0),
        head_code   = counter,
        body_code   = counter,
        tail_code   = counter,
        return_code = :(return count)
    ) |> eval
    
    @info "\tFor a pair of LongSequence{<:NucleicAcidAlphabet}s"
    
    counter = :(count += ambiguous_bitcount(x, y, A()))
    
    compile_2seq_bitpar(
        :count_ambiguous_bitpar,
        arguments = (:(seqa::LongSequence{A}), :(seqb::LongSequence{A})),
        parameters = (:(A<:NucleicAcidAlphabet),),
        init_code = :(count = 0),
        head_code = counter,
        body_code = counter,
        tail_code = counter,
        return_code = :(return count)
    ) |> eval
end


## For a single sequence.
# You can never have ambiguous bases in a 2-bit encoded nucleotide sequence.
Base.count(::typeof(isambiguous), seq::LongSequence{<:NucleicAcidAlphabet{2}}) = 0
Base.count(::typeof(isambiguous), seq::LongSequence{<:NucleicAcidAlphabet{4}}) = count_ambiguous_bitpar(seq)

## For a pair of sequences.
# A pair of 2-bit encoded sequences will never have ambiguous bases.
Base.count(::typeof(isambiguous), seqa::LongSequence{A}, seqb::LongSequence{A}) where {A<:NucleicAcidAlphabet{2}} = 0
Base.count(::typeof(isambiguous), seqa::LongSequence{A}, seqb::LongSequence{A}) where {A<:NucleicAcidAlphabet{4}} = count_ambiguous_bitpar(seqa, seqb)
Base.count(::typeof(isambiguous), seqa::LongSequence{<:NucleicAcidAlphabet{4}}, seqb::LongSequence{<:NucleicAcidAlphabet{2}}) = count(isambiguous, promote(seqa, seqb)...)
Base.count(::typeof(isambiguous), seqa::LongSequence{<:NucleicAcidAlphabet{2}}, seqb::LongSequence{<:NucleicAcidAlphabet{4}}) = count(isambiguous, promote(seqa, seqb)...)

# Counting certain sites
let
    @info "Compiling bit-parallel certainty counter for LongSequence{<:NucleicAcidAlphabet}"
    
    counter = :(count += certain_bitcount(x, y, A()))
    
    compile_2seq_bitpar(
        :count_certain_bitpar,
        arguments = (:(seqa::LongSequence{A}), :(seqb::LongSequence{A})),
        parameters = (:(A<:NucleicAcidAlphabet),),
        init_code = :(count = 0),
        head_code = counter,
        body_code = counter,
        tail_code = counter,
        return_code = :(return count)
    ) |> eval
end
Base.count(::typeof(iscertain), seqa::LongSequence{A}, seqb::LongSequence{A}) where {A<:NucleicAcidAlphabet{4}} = count_certain_bitpar(seqa, seqb)
Base.count(::typeof(iscertain), seqa::LongSequence{<:NucleicAcidAlphabet{4}}, seqb::LongSequence{<:NucleicAcidAlphabet{2}}) = count(iscertain, promote(seqa, seqb)...)
Base.count(::typeof(iscertain), seqa::LongSequence{<:NucleicAcidAlphabet{2}}, seqb::LongSequence{<:NucleicAcidAlphabet{4}}) = count(iscertain, promote(seqa, seqb)...)

# Counting gap sites
let
    @info "Compiling bit-parallel gap counter for LongSequence{<:NucleicAcidAlphabet}"
    
    counter = :(count += gap_bitcount(x, y, A()))
    
    count_empty = quote
        count += gap_bitcount(x, y, A())
        nempty = div(64, bits_per_symbol(A())) - div(offs, bits_per_symbol(A()))
        count -= nempty
    end
    
    compile_2seq_bitpar(
        :count_gap_bitpar,
        arguments = (:(seqa::LongSequence{A}), :(seqb::LongSequence{A})),
        parameters = (:(A<:NucleicAcidAlphabet),),
        init_code = :(count = 0),
        head_code = count_empty,
        body_code = counter,
        tail_code = count_empty,
        return_code = :(return count)
    ) |> eval
end
Base.count(::typeof(isgap), seqa::LongSequence{A}, seqb::LongSequence{A}) where {A<:NucleicAcidAlphabet{4}} = count_gap_bitpar(seqa, seqb)
Base.count(::typeof(isgap), seqa::LongSequence{<:NucleicAcidAlphabet{4}}, seqb::LongSequence{<:NucleicAcidAlphabet{2}}) = count(isgap, promote(seqa, seqb)...)
Base.count(::typeof(isgap), seqa::LongSequence{<:NucleicAcidAlphabet{2}}, seqb::LongSequence{<:NucleicAcidAlphabet{4}}) = count(isgap, promote(seqa, seqb)...)
