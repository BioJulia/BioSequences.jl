###
### String Decorators
###
###
### String literals for LongSequences
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

remove_newlines(s) = replace(s, r"\r|\n" => "")

"""
    @dna_str(seq, flag="s") -> LongDNA{4}

Create a `LongDNA{4}` sequence at parse time from string `seq`.
If `flag` is `"s"` ('static', the default), the sequence is created at parse time,
and inserted directly into the returned expression.
A static string ought not to be mutated
Alternatively, if `flag` is `"d"` (dynamic), a new sequence is parsed and created
whenever the code where is macro is placed is run.

See also: [`@aa_str`](@ref), [`@rna_str`](@ref)

# Examples
In the example below, the static sequence is created once, at parse time, NOT
when the function `f` is run. This means it is the _same_  sequence that is
pushed to repeatedly.
```jldoctest
julia> f() = dna"TAG";

julia> string(push!(f(), DNA_A)) # NB: Mutates static string!
"TAGA"

julia> string(push!(f(), DNA_A))
"TAGAA"

julia> f() = dna"TAG"d; # dynamically make seq

julia> string(push!(f(), DNA_A))
"TAGA"

julia> string(push!(f(), DNA_A))
"TAGA"
```
"""
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

"""
The `LongRNA{4}` equivalent to `@dna_str`

See also: [`@dna_str`](@ref), [`@aa_str`](@ref)

# Examples
```jldoctest
julia> rna"UCGUGAUGC"
9nt RNA Sequence:
UCGUGAUGC
```
"""
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


"""
The `AminoAcidAlphabet` equivalent to `@dna_str`

See also: [`@dna_str`](@ref), [`@rna_str`](@ref)

# Examples
```jldoctest
julia> aa"PKLEQC"
6aa Amino Acid Sequence:
PKLEQC
```
"""
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
