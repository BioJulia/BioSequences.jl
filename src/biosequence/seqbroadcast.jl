###
### Broadcasting with BioSequences.
###
### Customising the broadcasting machinery for BioSequences.
###
### This file is a part of BioJulia.
### License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md


# First we ensure that BioSequence types are actually passed to the broadcast
# machinery.
# 
# By default, this would otherwise `collect` the sequence into a `Vector` of
# `DNA` or `RNA` etc. This wastes time and memory, when broadcasting should
# just be able to use the BioSequence itself as if it were an array with a
# a single dimension, akin to Vector.
Base.broadcastable(x::BioSequence) = x



#=

Let's have a look at the example of broadcasting addition with a vector.

@less broadcast(+, [1,2,3], 1)
===
broadcast(f::Tf, As...) where {Tf} = materialize(broadcasted(f, As...))


@less Base.broadcasted(+, [1,2,3], 1)
===
@inline function broadcasted(f, arg1, arg2, args...)
    arg1′ = broadcastable(arg1)
    arg2′ = broadcastable(arg2)
    args′ = map(broadcastable, args)
    broadcasted(combine_styles(arg1′, arg2′, args′...), f, arg1′, arg2′, args′...)
end
@inline broadcasted(::S, f, args...) where S<:BroadcastStyle = Broadcasted{S}(f, args)

bc =  Base.broadcasted(Base.Broadcast.combine_styles(arg1, arg2), +, arg1, arg2)
Base.Broadcast.Broadcasted(+, ([1, 2, 3], 1))

typeof(bc)
Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Nothing, typeof(+), Tuple{Vector{Int64}, Int64}}


@less Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}}(+, (arg1, arg2))
===
function Broadcasted{Style}(f::F, args::Args, axes=nothing) where {Style, F, Args<:Tuple}
    # using Core.Typeof rather than F preserves inferrability when f is a type
    Broadcasted{Style, typeof(axes), Core.Typeof(f), Args}(f, args, axes)
end


@less materialize(bc)
===
@inline materialize(bc::Broadcasted) = copy(instantiate(bc))
materialize(x) = x


@less Base.Broadcast.instantiate(bc)
===
@inline function instantiate(bc::Broadcasted{Style}) where {Style}
    if bc.axes isa Nothing # Not done via dispatch to make it easier to extend instantiate(::Broadcasted{Style})
        axes = combine_axes(bc.args...)
    else
        axes = bc.axes
        check_broadcast_axes(axes, bc.args...)
    end
    return Broadcasted{Style}(bc.f, bc.args, axes)
end


@less Base.Broadcast.combine_axes(bc.args...)
===
@inline combine_axes(A, B...) = broadcast_shape(axes(A), combine_axes(B...))
@inline combine_axes(A, B) = broadcast_shape(axes(A), axes(B))


@less @less Base.Broadcast.broadcast_shape(axes([1,2,3]), axes(1))
===
broadcast_shape(shape::Tuple) = shape
broadcast_shape(shape::Tuple, shape1::Tuple, shapes::Tuple...) = broadcast_shape(_bcs(shape, shape1), shapes...)


@less Base.Broadcast._bcs(axes([1,2,3]), axes(1))
===
_bcs(shape::Tuple, ::Tuple{}) = (shape[1], _bcs(tail(shape), ())...)

_bcs consolidates two shapes into a single output shape

shape = axes([1,2,3])
@less Base.Broadcast._bcs(Base.tail(shape), ())
===
_bcs(::Tuple{}, ::Tuple{}) = ()
()

(shape[1], _bcs(tail(shape), ())...)


@less copy(Base.Broadcast.instantiate(bc))
===
@inline function copy(bc::Broadcasted{Style}) where {Style}
    ElType = combine_eltypes(bc.f, bc.args)
    if Base.isconcretetype(ElType)
        # We can trust it and defer to the simpler `copyto!`
        return copyto!(similar(bc, ElType), bc)
    end
    # When ElType is not concrete, use narrowing. Use the first output
    # value to determine the starting output eltype; copyto_nonleaf!
    # will widen `dest` as needed to accommodate later values.
    bc′ = preprocess(nothing, bc)
    iter = eachindex(bc′)
    y = iterate(iter)
    if y === nothing
        # if empty, take the ElType at face value
        return similar(bc′, ElType)
    end
    # Initialize using the first value
    I, state = y
    @inbounds val = bc′[I]
    dest = similar(bc′, typeof(val))
    @inbounds dest[I] = val
    # Now handle the remaining values
    # The typeassert gives inference a helping hand on the element type and dimensionality
    # (work-around for #28382)
    ElType′ = ElType === Union{} ? Any : ElType <: Type ? Type : ElType
    RT = dest isa AbstractArray ? AbstractArray{<:ElType′, ndims(dest)} : Any
    return copyto_nonleaf!(dest, bc′, iter, state, 1)::RT
end


@less Base.Broadcast.combine_eltypes(bc.f, bc.args)
===
combine_eltypes(f, args::Tuple) = promote_typejoin_union(Base._return_type(f, eltypes(args)))

@less Base.Broadcast.eltypes(bc.args)

@less Base._return_type(+, Base.Broadcast.eltypes(bc.args))

@less Base.Broadcast.promote_typejoin_union(Base._return_type(+, Base.Broadcast.eltypes(bc.args)))
===
function promote_typejoin_union(::Type{T}) where T
    if T === Union{}
        return Union{}
    elseif T isa UnionAll
        return Any # TODO: compute more precise bounds
    elseif T isa Union
        return promote_typejoin(promote_typejoin_union(T.a), promote_typejoin_union(T.b))
    elseif T <: Tuple
        return typejoin_union_tuple(T)
    else
        return T
    end
end


@less similar(bc, Base.Broadcast.combine_eltypes(bc.f, bc.args))
===
Base.similar(bc::Broadcasted, ::Type{T}) where {T} = similar(bc, T, axes(bc))
Base.similar(::Broadcasted{DefaultArrayStyle{N}}, ::Type{ElType}, dims) where {N,ElType} = similar(Array{ElType}, dims)
Base.similar(::Broadcasted{DefaultArrayStyle{N}}, ::Type{Bool}, dims) where N = similar(BitArray, dims)

@less copyto!(similar(bc, Base.Broadcast.combine_eltypes(bc.f, bc.args)), bc)
===
## general `copyto!` methods
# The most general method falls back to a method that replaces Style->Nothing
# This permits specialization on typeof(dest) without introducing ambiguities
@inline copyto!(dest::AbstractArray, bc::Broadcasted) = copyto!(dest, convert(Broadcasted{Nothing}, bc))


@less copyto!(dest, convert(Base.Broadcast.Broadcasted{Nothing}, bc))
===
@inline function copyto!(dest::AbstractArray, bc::Broadcasted{Nothing})
    axes(dest) == axes(bc) || throwdm(axes(dest), axes(bc))
    # Performance optimization: broadcast!(identity, dest, A) is equivalent to copyto!(dest, A) if indices match
    if bc.f === identity && bc.args isa Tuple{AbstractArray} # only a single input argument to broadcast!
        A = bc.args[1]
        if axes(dest) == axes(A)
            return copyto!(dest, A)
        end
    end
    bc′ = preprocess(dest, bc)
    # Performance may vary depending on whether `@inbounds` is placed outside the
    # for loop or not. (cf. https://github.com/JuliaLang/julia/issues/38086)
    @inbounds @simd for I in eachindex(bc′)
        dest[I] = bc′[I]
    end
    return dest
end


@less Base.Broadcast.preprocess(dest, bc)
===
# Preprocessing a `Broadcasted` does two things:
# * unaliases any arguments from `dest`
# * "extrudes" the arguments where it is advantageous to pre-compute the broadcasted indices
@inline preprocess(dest, bc::Broadcasted{Style}) where {Style} = Broadcasted{Style}(bc.f, preprocess_args(dest, bc.args), bc.axes)
preprocess(dest, x) = extrude(broadcast_unalias(dest, x))

@inline preprocess_args(dest, args::Tuple) = (preprocess(dest, args[1]), preprocess_args(dest, tail(args))...)
preprocess_args(dest, args::Tuple{Any}) = (preprocess(dest, args[1]),)
preprocess_args(dest, args::Tuple{}) = ()


Base.Broadcast.preprocess(dest, bc.args[1])
@less Base.Broadcast.broadcast_unalias(dest, bc.args[1])
===
# For broadcasted assignments like `broadcast!(f, A, ..., A, ...)`, where `A`
# appears on both the LHS and the RHS of the `.=`, then we know we're only
# going to make one pass through the array, and even though `A` is aliasing
# against itself, the mutations won't affect the result as the indices on the
# LHS and RHS will always match. This is not true in general, but with the `.op=`
# syntax it's fairly common for an argument to be `===` a source.
broadcast_unalias(dest, src) = dest === src ? src : unalias(dest, src)
broadcast_unalias(::Nothing, src) = src


@less Base.Broadcast.unalias(dest, bc.args[1])
===
## rudimentary aliasing detection ##
"""
    Base.unalias(dest, A)

Return either `A` or a copy of `A` in a rough effort to prevent modifications to `dest` from
affecting the returned object. No guarantees are provided.

Custom arrays that wrap or use fields containing arrays that might alias against other
external objects should provide a [`Base.dataids`](@ref) implementation.

This function must return an object of exactly the same type as `A` for performance and type
stability. Mutable custom arrays for which [`copy(A)`](@ref) is not `typeof(A)` should
provide a [`Base.unaliascopy`](@ref) implementation.

See also [`Base.mightalias`](@ref).
"""
unalias(dest, A::AbstractArray) = mightalias(dest, A) ? unaliascopy(A) : A
unalias(dest, A::AbstractRange) = A
unalias(dest, A) = A

"""
    Base.mightalias(A::AbstractArray, B::AbstractArray)

Perform a conservative test to check if arrays `A` and `B` might share the same memory.

By default, this simply checks if either of the arrays reference the same memory
regions, as identified by their [`Base.dataids`](@ref).
"""
mightalias(A::AbstractArray, B::AbstractArray) = !isbits(A) && !isbits(B) && !_isdisjoint(dataids(A)
, dataids(B))
mightalias(x, y) = false

"""
    Base.dataids(A::AbstractArray)

Return a tuple of `UInt`s that represent the mutable data segments of an array.

Custom arrays that would like to opt-in to aliasing detection of their component
parts can specialize this method to return the concatenation of the `dataids` of
their component parts.  A typical definition for an array that wraps a parent is
`Base.dataids(C::CustomArray) = dataids(C.parent)`.
"""
dataids(A::AbstractArray) = (UInt(objectid(A)),)
dataids(A::Array) = (UInt(pointer(A)),)
dataids(::AbstractRange) = ()
dataids(x) = ()

=#

