
@inline function inbounds_getindex(x::ShortSequence, i::Integer)
    return reinterpret(eltype(x), 0x01 << ((encoded_data(x) >> 2(length(x) - i)) & 0b11))
end