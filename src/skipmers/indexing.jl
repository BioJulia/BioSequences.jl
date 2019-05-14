
@inline offset(::Type{S}, i::Integer) where {S <: Skipmer} = 2(kmersize(S) - i)
@inline offset(x::Skipmer, i::Integer) = offset(typeof(x), i)

@inline function inbounds_getindex(x::Skipmer, i::Integer)
    return reinterpret(eltype(x), 0x01 << ((encoded_data(x) >> offset(x, i)) & 0b11))
end

