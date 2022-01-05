```@meta
CurrentModule = BioSequences
DocTestSetup = quote
    using BioSequences
end
```

# Iteration

As you might expect, sequence types are iterators over their elements:

```jldoctest
julia> n = 0
0

julia> for nt in dna"ATNGNNT"
           if nt == DNA_N
               global n += 1
           end
       end

julia> n
3

```