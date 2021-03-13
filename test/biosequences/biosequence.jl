# Create new biosequence for testing
struct SimpleSeq <: BioSequence{RNAAlphabet{2}}
    x::Vector{UInt}
end

SimpleSeq(it) = SimpleSeq([BioSequences.encode(Alphabet(SimpleSeq), i) for i in it])
Base.copy(x::SimpleSeq) = SimpleSeq(copy(x.x))
Base.length(x::SimpleSeq) = length(x.x)
BioSequences.encoded_data_eltype(::Type{SimpleSeq}) = UInt
BioSequences.extract_encoded_element(x::SimpleSeq, i::Integer) = x.x[i]

BioSequences.encoded_setindex!(x::SimpleSeq, e::UInt, i::Integer) = x.x[i] = e
SimpleSeq(::UndefInitializer, x::Integer) = SimpleSeq(zeros(UInt, x))
resize!(x::SimpleSeq, len::Int) = (resize!(x.x, len), x)

@testset "Basics" begin
    seq = SimpleSeq([RNA_C, RNA_G, RNA_U])
    
    @test seq isa BioSequence{RNAAlphabet{2}}
    @test Alphabet(seq) === RNAAlphabet{2}()
    @test empty(seq) isa BioSequence
    @test eltype(typeof(seq)) == RNA
    @test eltype(seq) == RNA

    @test copy(seq) == seq
    @test copy(seq) !== seq
    
    @test length(seq) == 3
    @test size(seq) == (3,)
    @test length(empty(seq)) == 0
    @test isempty(empty(seq))
    @test isempty(empty(seq))

    @test firstindex(seq) == 1
    @test lastindex(seq) == 3
    @test collect(eachindex(seq)) == [1, 2, 3]
    @test collect(keys(seq)) == [1, 2, 3]

    @test nextind(seq, firstindex(seq)) == 2
    @test prevind(seq, lastindex(seq)) == 2
    @test prevind(seq, 2) == 1

    @test_throws EncodeError SimpleSeq([RNA_C, RNA_G, RNA_M])
    @test_throws EncodeError SimpleSeq([RNA_Gap])
    @test_throws MethodError SimpleSeq([DNA_C, DNA_G, DNA_C])
    @test_throws MethodError SimpleSeq(1:3)
end
