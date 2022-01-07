###
### LongSequence specific specializations of src/biosequence/counting.jl
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

# Counting GC positions
let
    counter = :(n += gc_bitcount(chunk, Alphabet(seq)))
    compile_bitpar(
        :count_gc_bitpar,
        arguments   = (:(seq::SeqOrView{<:NucleicAcidAlphabet}),),
        init_code   = :(n = 0),
        head_code   = counter,
        body_code   = counter,
        tail_code   = counter,
        return_code = :(return n % Int)
    ) |> eval
end

Base.count(::typeof(isGC), seq::SeqOrView{<:NucleicAcidAlphabet}) = count_gc_bitpar(seq)

# Counting mismatches
let
    counter = :(count += mismatch_bitcount(x, y, A()))
    
    compile_2seq_bitpar(
        :count_mismatches_bitpar,
        arguments = (:(seqa::SeqOrView{A}), :(seqb::SeqOrView{A})),
        parameters = (:(A<:NucleicAcidAlphabet),),
        init_code = :(count = 0),
        head_code = counter,
        body_code = counter,
        tail_code = counter,
        return_code = :(return count % Int)
    ) |> eval
end
Base.count(::typeof(!=), seqa::SeqOrView{A}, seqb::SeqOrView{A}) where {A<:NucleicAcidAlphabet} = count_mismatches_bitpar(seqa, seqb)
Base.count(::typeof(!=), seqa::NucleicSeqOrView, seqb::NucleicSeqOrView) = count(!=, promote(seqa, seqb)...)

# Counting matches
let
    counter = :(count += match_bitcount(x, y, A()))
    
    count_empty = quote
        count += match_bitcount(x, y, A())
        nempty = div(64, bits_per_symbol(A())) - div(Int(offs), bits_per_symbol(A()))
        count -= nempty
    end
    
    compile_2seq_bitpar(
        :count_matches_bitpar,
        arguments = (:(seqa::SeqOrView{A}), :(seqb::SeqOrView{A})),
        parameters = (:(A<:NucleicAcidAlphabet),),
        init_code = :(count = 0),
        head_code = count_empty,
        body_code = counter,
        tail_code = count_empty,
        return_code = :(return count % Int)
    ) |> eval
end
Base.count(::typeof(==), seqa::SeqOrView{A}, seqb::SeqOrView{A}) where {A<:NucleicAcidAlphabet} = count_matches_bitpar(seqa, seqb)
Base.count(::typeof(==), seqa::NucleicSeqOrView, seqb::NucleicSeqOrView) = count(==, promote(seqa, seqb)...)

# Counting ambiguous sites
# ------------------------
let
    counter = :(count += ambiguous_bitcount(chunk, Alphabet(seq)))
    
    compile_bitpar(
        :count_ambiguous_bitpar,
        arguments   = (:(seq::SeqOrView{<:NucleicAcidAlphabet}),),
        init_code   = :(count = 0),
        head_code   = counter,
        body_code   = counter,
        tail_code   = counter,
        return_code = :(return count % Int)
    ) |> eval
    
    counter = :(count += ambiguous_bitcount(x, y, A()))
    
    compile_2seq_bitpar(
        :count_ambiguous_bitpar,
        arguments = (:(seqa::SeqOrView{A}), :(seqb::SeqOrView{A})),
        parameters = (:(A<:NucleicAcidAlphabet),),
        init_code = :(count = 0),
        head_code = counter,
        body_code = counter,
        tail_code = counter,
        return_code = :(return count % Int)
    ) |> eval
end


## For a single sequence.
# You can never have ambiguous bases in a 2-bit encoded nucleotide sequence.
Base.count(::typeof(isambiguous), seq::SeqOrView{<:NucleicAcidAlphabet{2}}) = 0
Base.count(::typeof(isambiguous), seq::SeqOrView{<:NucleicAcidAlphabet{4}}) = count_ambiguous_bitpar(seq)

## For a pair of sequences.
# A pair of 2-bit encoded sequences will never have ambiguous bases.
Base.count(::typeof(isambiguous), seqa::SeqOrView{A}, seqb::SeqOrView{A}) where {A<:NucleicAcidAlphabet{2}} = 0
Base.count(::typeof(isambiguous), seqa::SeqOrView{A}, seqb::SeqOrView{A}) where {A<:NucleicAcidAlphabet{4}} = count_ambiguous_bitpar(seqa, seqb)
Base.count(::typeof(isambiguous), seqa::SeqOrView{<:NucleicAcidAlphabet{4}}, seqb::SeqOrView{<:NucleicAcidAlphabet{2}}) = count(isambiguous, promote(seqa, seqb)...)
Base.count(::typeof(isambiguous), seqa::SeqOrView{<:NucleicAcidAlphabet{2}}, seqb::SeqOrView{<:NucleicAcidAlphabet{4}}) = count(isambiguous, promote(seqa, seqb)...)

# Counting certain sites
let
    counter = :(count += certain_bitcount(x, y, A()))
    
    compile_2seq_bitpar(
        :count_certain_bitpar,
        arguments = (:(seqa::SeqOrView{A}), :(seqb::SeqOrView{A})),
        parameters = (:(A<:NucleicAcidAlphabet),),
        init_code = :(count = 0),
        head_code = counter,
        body_code = counter,
        tail_code = counter,
        return_code = :(return count % Int)
    ) |> eval
end
Base.count(::typeof(iscertain), seqa::SeqOrView{A}, seqb::SeqOrView{A}) where {A<:NucleicAcidAlphabet{4}} = count_certain_bitpar(seqa, seqb)
Base.count(::typeof(iscertain), seqa::SeqOrView{<:NucleicAcidAlphabet{4}}, seqb::SeqOrView{<:NucleicAcidAlphabet{2}}) = count(iscertain, promote(seqa, seqb)...)
Base.count(::typeof(iscertain), seqa::SeqOrView{<:NucleicAcidAlphabet{2}}, seqb::SeqOrView{<:NucleicAcidAlphabet{4}}) = count(iscertain, promote(seqa, seqb)...)

# Counting gap sites
let
    count_empty = quote
        Alph = Alphabet(seq)
        count += gap_bitcount(chunk, Alph)
        count -= div(n_bits_masked, bits_per_symbol(Alph))
    end
    counter = :(count += gap_bitcount(chunk, Alphabet(seq)))
    
    compile_bitpar(
        :count_gap_bitpar,
        arguments   = (:(seq::SeqOrView{<:NucleicAcidAlphabet{4}}),),
        init_code   = :(count = 0),
        head_code   = count_empty,
        body_code   = counter,
        tail_code   = count_empty,
        return_code = :(return count % Int)
    ) |> eval

    count_empty = quote
        Alph = Alphabet(seqa)
        count += gap_bitcount(x, y, Alph)
        nempty = div(64, bits_per_symbol(Alph)) - div(offs, bits_per_symbol(Alph))
        count -= nempty
    end
    counter = :(count += gap_bitcount(x, y, A()))
    
    compile_2seq_bitpar(
        :count_gap_bitpar,
        arguments = (:(seqa::SeqOrView{A}), :(seqb::SeqOrView{A})),
        parameters = (:(A<:NucleicAcidAlphabet),),
        init_code = :(count = 0),
        head_code = count_empty,
        body_code = counter,
        tail_code = count_empty,
        return_code = :(return count % Int)
    ) |> eval
end
Base.count(::typeof(isgap), seqa::SeqOrView{A}, seqb::SeqOrView{A}) where {A<:NucleicAcidAlphabet{4}} = count_gap_bitpar(seqa, seqb)
Base.count(::typeof(isgap), seqa::SeqOrView{A}) where {A<:NucleicAcidAlphabet{4}} = count_gap_bitpar(seqa)
Base.count(::typeof(isgap), seqa::SeqOrView{<:NucleicAcidAlphabet{4}}, seqb::SeqOrView{<:NucleicAcidAlphabet{2}}) = count(isgap, promote(seqa, seqb)...)
Base.count(::typeof(isgap), seqa::SeqOrView{<:NucleicAcidAlphabet{2}}, seqb::SeqOrView{<:NucleicAcidAlphabet{4}}) = count(isgap, promote(seqa, seqb)...)
