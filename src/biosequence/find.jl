###
### Finders
###

function Base.findnext(val, seq::BioSequence, start::Integer)
    checkbounds(seq, start)
    v = convert(eltype(seq), val)
    for i in start:lastindex(seq)
        if inbounds_getindex(seq, i) == v
            return i
        end
    end
    return nothing
end

function Base.findnext(f::Function, seq::BioSequence, start::Integer)
    checkbounds(seq, start)
    for i in start:lastindex(seq)
        if f(inbounds_getindex(seq, i))
            return i
        end
    end
    return nothing
end

# No ambiguous sites can exist in a nucleic acid sequence using the two-bit alphabet.
Base.findnext(::typeof(isambiguous), seq::BioSequence{<:NucleicAcidAlphabet{2}}, from::Integer) = nothing

function Base.findprev(val, seq::BioSequence, start::Integer)
    checkbounds(seq, start)
    v = convert(eltype(seq), val)
    for i in start:-1:1
        if inbounds_getindex(seq, i) == v
            return i
        end
    end
    return nothing
end

function Base.findprev(f::Function, seq::BioSequence, start::Integer)
    checkbounds(seq, start)
    for i in start:-1:firstindex(seq)
        if f(inbounds_getindex(seq, i))
            return i
        end
    end
    return nothing
end

Base.findfirst(val, seq::BioSequence) = findnext(val, seq, 1)
Base.findlast(val, seq::BioSequence) = findprev(val, seq, lastindex(seq))

