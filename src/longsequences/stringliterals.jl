###
### String Decorators
###
###
### String literals for LongSequences
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

remove_newlines(s) = replace(s, r"\r|\n" => "")

macro dna_str(seq, flag)
    if flag == "s"
        return LongDNA{4}(remove_newlines(seq))
    elseif flag == "d"
        return quote
            LongDNA{4}($(remove_newlines(seq)))
        end
    end
    error("Invalid DNA flag: '$(flag)'")
end

macro dna_str(seq)
    return LongDNA{4}(remove_newlines(seq))
end

macro rna_str(seq, flag)
    if flag == "s"
        return LongRNA{4}(remove_newlines(seq))
    elseif flag == "d"
        return quote
            LongRNA{4}($(remove_newlines(seq)))
        end
    end
    error("Invalid RNA flag: '$(flag)'")
end

macro rna_str(seq)
    return LongRNA{4}(remove_newlines(seq))
end

macro aa_str(seq, flag)
    if flag == "s"
        return LongAA(remove_newlines(seq))
    elseif flag == "d"
        return quote
            LongAA($(remove_newlines(seq)))
        end
    end
    error("Invalid Amino Acid flag: '$(flag)'")
end

macro aa_str(seq)
    return LongAA(remove_newlines(seq))
end
