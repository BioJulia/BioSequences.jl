@inline _offset(K, i) = 2(K - i)
@inline offset(::Type{T}, i::Integer) where {A,K,T<:AbstractMer{A,K}} = 2(K - i)
@inline offset(x::AbstractMer, i::Integer) = offset(typeof(x), i)

@inline function inbounds_getindex(x::AbstractMer, i::Integer)
    return reinterpret(eltype(x), 0x01 << ((encoded_data(x) >> offset(x, i)) & 0b11))
end

