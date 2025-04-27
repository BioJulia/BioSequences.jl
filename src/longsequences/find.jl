# Zeros out all nibbles except the ones with 4bit A,C,G or T
@inline function iscertain_kernel(::NucleicAcidAlphabet{4}, x::UInt64)
    x = enumerate_nibbles(x)
    y = (x & 0x4444444444444444) >> 2
    y |= (x & 0x2222222222222222) >> 1
    return x & ~y & 0x1111111111111111
end

# Zeros out all which is not a certain (normal or AA_Term)
@inline function iscertain_kernel(::AminoAcidAlphabet, x::UInt64)
    # 1. Set normal to FF, others to 00
    y = simd_lt_byte(x, 0x16)

    # 2. Set Term to FF, others to 00
    z = set_zero_encoding(BitsPerSymbol{8}(), x ⊻ 0x1a1a1a1a1a1a1a1a) * 0xFF

    # 3: OR them
    return y | z
end

function _findnext(::typeof(iscertain), seq::SeqOrView{<:Union{<:NucleicAcidAlphabet{4}, AminoAcidAlphabet}}, i::Int)
    f = @inline u -> iscertain_kernel(Alphabet(seq), u)
    return _findnext_nonzero(f, seq, i)
end

function _findprev(::typeof(iscertain), seq::SeqOrView{<:Union{<:NucleicAcidAlphabet{4}, AminoAcidAlphabet}}, i::Int)
    f = @inline u -> iscertain_kernel(Alphabet(seq), u)
    return _findprev_nonzero(f, seq, i)
end

@inline function simd_lt_byte(x::UInt64, byte::UInt8)
    T = NTuple{8, VecElement{UInt8}}
    x = reinterpret(T, x)
    y = ntuple(i -> VecElement(byte), Val{8}())
    s = """
    %res = icmp ult <8 x i8> %0, %1
    %resb = sext <8 x i1> %res to <8 x i8>
    ret <8 x i8> %resb
    """
    z = Core.Intrinsics.llvmcall(s, NTuple{8, VecElement{UInt8}}, Tuple{NTuple{8, VecElement{UInt8}}, NTuple{8, VecElement{UInt8}}}, x, y)
    return reinterpret(UInt64, z)
end

@inline function sub_byte(x::UInt64, byte::UInt8)
    T = NTuple{8, VecElement{UInt8}}
    x = reinterpret(T, x)
    y = ntuple(i -> VecElement(byte), Val{8}())
    s = """
    %res = sub <8 x i8> %0, %1
    ret <8 x i8> %res
    """
    z = Core.Intrinsics.llvmcall(s, NTuple{8, VecElement{UInt8}}, Tuple{NTuple{8, VecElement{UInt8}}, NTuple{8, VecElement{UInt8}}}, x, y)
    return reinterpret(UInt64, z)
end

# Zeros out all nibbles encoding 4bit A,C,G or T
@inline function uncertain_kernel(::NucleicAcidAlphabet{4}, x::UInt64)
    return enumerate_nibbles(x) ⊻ 0x1111111111111111
end

# Zero out normal AAs and AA_Term
@inline function uncertain_kernel(::AminoAcidAlphabet, x::UInt64)
    # Zero out normal AA, set rest to 0xFF
    y = ~simd_lt_byte(x, 0x16)

    # Zero out 0x1a, set rest to other bitpatterns
    z = x ⊻ 0x1a1a1a1a1a1a1a1a

    return y & z
end

function _findnext(
        ::ComposedFunction{typeof(!), typeof(iscertain)},
        seq::SeqOrView{<:Union{<:NucleicAcidAlphabet{4}, AminoAcidAlphabet}},
        i::Int,
    )
    f = @inline u -> uncertain_kernel(Alphabet(seq), u)
    return _findnext_nonzero(f, seq, i)
end

function _findprev(
        ::ComposedFunction{typeof(!), typeof(iscertain)},
        seq::SeqOrView{<:Union{<:NucleicAcidAlphabet{4}, AminoAcidAlphabet}},
        i::Int,
    )
    f = @inline u -> uncertain_kernel(Alphabet(seq), u)
    return _findprev_nonzero(f, seq, i)
end

# Zeros out A, C, G, T or Gap
@inline function ambiguous_kernel(A::NucleicAcidAlphabet{4}, x::UInt64)
    # The y part makes every nibble 0xF, unless it's 0 to begin with
    y = x | (x >>> 2)
    y |= y >>> 1
    y &= 0x1111111111111111
    y *= 0x0F
    return uncertain_kernel(A, x) & y
end

# Zero out all except ambiguous symbols (AA_B, AA_J, AA_Z, AA_X), 0x16:0x19
@inline function ambiguous_kernel(::AminoAcidAlphabet, x::UInt64)
    return simd_lt_byte(sub_byte(x, 0x16), 0x04)
end

function _findnext(::typeof(isambiguous), seq::SeqOrView{<:Union{<:NucleicAcidAlphabet{4}, AminoAcidAlphabet}}, i::Int)
    f = @inline u -> ambiguous_kernel(Alphabet(seq), u)
    return _findnext_nonzero(f, seq, i)
end

function _findprev(::typeof(isambiguous), seq::SeqOrView{<:Union{<:NucleicAcidAlphabet{4}, AminoAcidAlphabet}}, i::Int)
    f = @inline u -> ambiguous_kernel(Alphabet(seq), u)
    return _findprev_nonzero(f, seq, i)
end

# Zeros out all except A, C, G, T and Gap
@inline function unambiguous_kernel(::NucleicAcidAlphabet{4}, x::UInt64)
    y = enumerate_nibbles(x)
    y = (y >>> 1) | (y >>> 2)
    y &= 0x1111111111111111
    return y ⊻= 0x1111111111111111
end

# Zero out the four ambiguous amino acids B, J, Z, X
@inline function unambiguous_kernel(::AminoAcidAlphabet, x::UInt64)
    return ~simd_lt_byte(sub_byte(x, 0x16), 0x04)
end

function _findnext(
        ::ComposedFunction{typeof(!), typeof(isambiguous)},
        seq::SeqOrView{<:Union{<:NucleicAcidAlphabet{4}, AminoAcidAlphabet}},
        i::Int,
    )
    f = @inline u -> unambiguous_kernel(Alphabet(seq), u)
    return _findnext_nonzero(f, seq, i)
end

function _findprev(
        ::ComposedFunction{typeof(!), typeof(isambiguous)},
        seq::SeqOrView{<:Union{<:NucleicAcidAlphabet{4}, AminoAcidAlphabet}},
        i::Int,
    )
    f = @inline u -> unambiguous_kernel(Alphabet(seq), u)
    return _findprev_nonzero(f, seq, i)
end

# For debug
#=
function make_all_nibbles()
    x = UInt64(0)
    for i in 0:15
        x = (x << 4) | i
    end
    x
end
=#
