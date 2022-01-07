module TestBioSequences

using Test
using Random
using StableRNGs
using LinearAlgebra: normalize
import BioSymbols
using BioSequences
using StatsBase
using YAML

# Test for the absence of ubound type parameters in the package
@test length(Test.detect_unbound_args(BioSequences, recursive=true)) == 0

# Test utils not dependent on BioSymbols
include("utils.jl")

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
