module TestExt

using Test: Test, GenericString
using BioSequences: BioSequence

# Resolve ambiguity with BioSequences' generic AbstractString constructor.
# Can be removed when that constructor is removed in the next breaking release.
Test.GenericString(seq::BioSequence) = GenericString(String(seq))

end # module
