module TestBioSequences

using Test

using Random
using StableRNGs
using LinearAlgebra: normalize
using BioSymbols
using BioSequences
using StatsBase
using YAML
using Aqua

@testset "Aqua" begin
    # On some Julia versions (e.g. 1.10), Core.Compiler has its own LazyString, which
    # is ambiguous with BioSequences' generic AbstractString constructor.
    # This can't be resolved without pirating Core.Compiler.
    # Remove this exclusion when the generic constructor is removed.
    ambiguity_exclude = Any[]
    if isdefined(Core.Compiler, :LazyString) && Core.Compiler.LazyString !== LazyString
        push!(ambiguity_exclude, Core.Compiler.LazyString)
    end
    Aqua.test_all(BioSequences; ambiguities=(exclude=ambiguity_exclude,))
end

# Test for the absence of ubound type parameters in the package
@test length(Test.detect_unbound_args(BioSequences, recursive=true)) == 0

# Test utils not dependent on BioSymbols
include("utils.jl")
include("bit-manipulation.jl")

@testset "Alphabet" begin
    include("alphabet.jl")
end

@testset "BioSequences" begin
    include("biosequences/biosequence.jl")
    include("biosequences/indexing.jl")
    include("biosequences/misc.jl")
end

@testset "LongSequences" begin
    include("longsequences/basics.jl")
    include("longsequences/conversion.jl")
    include("longsequences/seqview.jl")
    include("longsequences/hashing.jl")
    include("longsequences/iteration.jl")
    include("longsequences/mutability.jl")
    include("longsequences/print.jl")
    include("longsequences/transformations.jl")
    include("longsequences/predicates.jl")
    include("longsequences/find.jl")
    include("longsequences/randseq.jl")
    include("longsequences/shuffle.jl")
end
include("translation.jl")
include("counting.jl")

@testset "Search" begin
    include("search/ExactSearchQuery.jl")
    include("search/ApproximateSearchQuery.jl")
    include("search/regex.jl")
    include("search/pwm.jl")
end

end
