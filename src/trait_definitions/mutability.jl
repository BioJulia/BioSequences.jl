
"""
# Mutable and Immutable `BioSequence` types.

Concrete subtypes of `BioSequence` can be mutable or immutable.
"""
abstract type Mutability end

struct MutableSequence end

struct ImmutableSequence end

function Mutability end
