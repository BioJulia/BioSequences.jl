
"""
Return the `Alpahbet` type defining the possible biological symbols 
and their encoding for a given biological sequence.
"""
function alphabet_t end

@inline function alphabet_t(seq::BioSequence)
    return alphabet_t(typeof(seq))
end
    


