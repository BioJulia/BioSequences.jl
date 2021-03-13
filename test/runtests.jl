module TestBioSequences

using Test
using Random
using StableRNGs
using LinearAlgebra: normalize
import BioSymbols
using BioSequences
using StatsBase
using YAML

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
    #=
    include("longsequences/hashing.jl")
    include("longsequences/iteration.jl")
    include("longsequences/subseq.jl")
    include("longsequences/mutability.jl")
    include("longsequences/print.jl")
    include("longsequences/transformations.jl")
    include("longsequences/mutability.jl")
    include("longsequences/predicates.jl")
    include("longsequences/find.jl")
    include("longsequences/randseq.jl")
    include("longsequences/shuffle.jl")
    =#
end

@testset "Search" begin
    include("search/exact.jl")
    include("search/approximate.jl")
    include("search/regex.jl")
    include("search/pwm.jl")
end

include("counting.jl")
include("translation.jl")
=#

end
