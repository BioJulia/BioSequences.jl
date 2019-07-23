###
### LongSequence specific specializations of src/biosequence/counting.jl
###

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