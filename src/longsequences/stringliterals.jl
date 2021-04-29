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
        return LongDNASeq(remove_newlines(seq))
    elseif flag == "d"
        return quote
            LongDNASeq($(remove_newlines(seq)))
        end
    end
    error("Invalid DNA flag: '$(flag)'")
end

macro dna_str(seq)
    return LongDNASeq(remove_newlines(seq))
end

macro rna_str(seq, flag)
    if flag == "s"
        return LongRNASeq(remove_newlines(seq))
    elseif flag == "d"
        return quote
            LongRNASeq($(remove_newlines(seq)))
        end
    end
    error("Invalid RNA flag: '$(flag)'")
end

macro rna_str(seq)
    return LongRNASeq(remove_newlines(seq))
end

macro aa_str(seq, flag)
    if flag == "s"
        return LongAminoAcidSeq(remove_newlines(seq))
    elseif flag == "d"
        return quote
            LongAASeq($(remove_newlines(seq)))
        end
    end
    error("Invalid Amino Acid flag: '$(flag)'")
end

macro aa_str(seq)
    return LongAminoAcidSeq(remove_newlines(seq))
end
