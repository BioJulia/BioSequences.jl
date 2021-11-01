###
### Mer specific specializations of src/biosequence/indexing.jl
###

@inline _offset(K, i) = 2(K - i)
@inline offset(::Type{T}, i::Integer) where {A,K,T<:AbstractMer{A,K}} = 2(K - i)
@inline offset(x::AbstractMer, i::Integer) = offset(typeof(x), i)

@inline function inbounds_getindex(x::AbstractMer, i::Integer)
    return reinterpret(eltype(x), 0x01 << ((encoded_data(x) >> offset(x, i)) & 0b11))
end

function Base.getindex(mer::T, r::OrdinalRange{<:Integer, <:Integer}) where {A,K,T<:AbstractMer{A,K}}
    parameterless_mer_type = Base.typename(T).wrapper
    new_K = length(r)
    new_T = parameterless_mer_type{A, new_K}
    return new_T(mer[i] for i in r)
end
