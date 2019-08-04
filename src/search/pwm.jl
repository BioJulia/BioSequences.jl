# Motif Search based on Position Weight Matrix
# ============================================
#
# This file defines two types of matrices: PFM and PWM.  PFM is a thin wrapper
# of Matrix and stores non-negative frequency values for each symbol and
# position. PWM is a position weighted matrix and stores score values like PFM.
# PFMs are mutable like usual matrices but PWMs cache partial scores to each
# position and hence must be immutable so that updating operations won't cause
# inconsistency.  PFM is created from a raw frequency matrix or a set of
# sequences. PWM is created from a raw score matrix or a PFM.
#
# The motif search algorithm implemented here uses PWM. It tries to find the
# nearest position from which the total score exceeds a specified threshold. The
# algorithm stops searching immediately when it finds that the total score will
# not exceeds the threshold even if the maximum partial score is achieved from
# that position. This uses the partial score cache stored in a PWM object.

"""
Position frequency matrix.
"""
struct PFM{S,T<:Real} <: AbstractMatrix{T}
    # row: symbol, column: position
    data::Matrix{T}

    function PFM{S,T}(m::AbstractMatrix{T}) where {S<:Union{DNA,RNA},T}
        if size(m, 1) != 4
            throw(ArgumentError("PFM must have four rows"))
        end
        return new{S,T}(m)
    end
end

function PFM{S}(m::AbstractMatrix{T}) where {S, T}
    return PFM{S,T}(m)
end
Base.convert(::Type{PFM{S}}, m::AbstractMatrix{T}) where {S, T} = PFM{S}(m)

function Base.convert(::Type{Matrix{T}}, m::PFM) where T
    return convert(Matrix{T}, m.data)
end

function PFM(set)
    return PFM(collect(set))
end

function PFM(set::Vector)
    if isempty(set)
        throw(ArgumentError("sequence set must be non-empty"))
    end
    S = eltype(eltype(set))
    if S ∉ (DNA, RNA)
        throw(ArgumentError("sequence element must be DNA or RNA"))
    end
    len = length(set[1])
    freq = zeros(Int, (4, len))
    for i in 1:lastindex(set)
        seq = set[i]
        if eltype(seq) != S
            throw(ArgumentError("sequence element must be $(S)"))
        elseif length(seq) != len
            throw(ArgumentError("all sequences must be of the same length"))
        end
        for j in 1:len
            s = seq[j]
            if iscertain(s)  # ignore uncertain symbols
                freq[index_nuc(s, j)] += 1
            end
        end
    end
    return PFM{S}(freq)
end

# Broadcasting
struct PFMBroadcastStyle{S} <: Broadcast.BroadcastStyle end
Base.BroadcastStyle(::Type{PFM{S,T}}) where {S,T} = PFMBroadcastStyle{S}()
Base.BroadcastStyle(s1::PFMBroadcastStyle, s2::Base.BroadcastStyle) where {S,T} = s1
function Base.similar(bc::Broadcast.Broadcasted{PFMBroadcastStyle{S}}, elt::Type{T}) where {S, T}
    return PFM{S, T}(similar(Array{T}, axes(bc)))
end

function Base.IndexStyle(::Type{<:PFM})
    return IndexLinear()
end

function Base.size(m::PFM)
    return size(m.data)
end

function Base.getindex(m::PFM, i::Integer)
    return m.data[i]
end

function Base.getindex(m::PFM{S}, s::S, j::Integer) where S<:Union{DNA,RNA}
    return m.data[index_nuc(s, j)]
end

function Base.setindex!(m::PFM, val, i::Integer)
    return setindex!(m.data, val, i)
end

function Base.setindex!(m::PFM, val, i::Integer, j::Integer)
    return setindex!(m.data, val, i, j)
end

function Base.setindex!(m::PFM, val, s::S, j::Integer) where S<:Union{DNA,RNA}
    return setindex!(m.data, val, index_nuc(s, j))
end

function Base.show(io::IO, m::PFM{<:Union{DNA,RNA}})
    show_nuc_matrix(io, m)
end

function Base.show(io::IO, ::MIME"text/plain", m::PFM{<:Union{DNA,RNA}})
    show_nuc_matrix(io, m)
end


"""
Position weight matrix.
"""
struct PWM{S,T<:Real} <: AbstractMatrix{T}
    # symbol type
    symtype::Type{S}

    # score matrix (row: symbols, column: position)
    data::Matrix{T}

    # maxscore[i] == sum(maximum(data, 1)[i:end])
    maxscore::Vector{T}

    function PWM{S,T}(pwm::AbstractMatrix) where {S<:Union{DNA,RNA},T}
        if size(pwm, 1) != 4
            throw(ArgumentError("PWM must have four rows"))
        end
        # make a copy for safety
        return new{S,T}(S, copy(pwm), make_maxscore(pwm))
    end
end

"""
    PWM(pfm::PFM{<:Union{DNA,RNA}}; prior=fill(1/4,4))

Create a position weight matrix from a position frequency matrix `pfm`.

The positive weight matrix will be `log2.((pfm ./ sum(pfm, 1)) ./ prior)`.
"""
function PWM(pfm::PFM{S}; prior=fill(1/4, 4)) where S <: Union{DNA,RNA}
    if !all(x -> x > 0, prior)
        throw(ArgumentError("prior must be positive"))
    elseif sum(prior) ≉ 1
        throw(ArgumentError("prior must be sum to 1"))
    end
    prob = pfm ./ sum(pfm, dims=1)
    return PWM{S,Float64}(log2.(prob ./ prior))
end

"""
    PWM{T}(pwm::AbstractMatrix)

Create a PWM from a raw matrix `pwm` (rows are symbols and columns are
positions).

Examples
--------
```julia-repl
julia> pwm = [
           0.065 -0.007 0.042 0.016 0.003  # A
           0.025 -0.007 0.016 0.058 0.082  # C
           0.066  0.085 0.070 0.054 0.010  # G
           0.016  0.025 0.051 0.058 0.038  # T
       ]
4×5 Array{Float64,2}:
 0.065  -0.007  0.042  0.016  0.003
 0.025  -0.007  0.016  0.058  0.082
 0.066   0.085  0.07   0.054  0.01
 0.016   0.025  0.051  0.058  0.038

julia> PWM{DNA}(pwm)
BioSequences.PWM{BioSymbols.DNA,Float64}:
 A  0.065 -0.007  0.042  0.016  0.003
 C  0.025 -0.007  0.016  0.058  0.082
 G  0.066  0.085  0.07   0.054  0.01
 T  0.016  0.025  0.051  0.058  0.038

```
"""
function PWM{T}(pwm::AbstractMatrix) where T <: Union{DNA,RNA}
    return PWM{T,eltype(pwm)}(pwm)
end

function Base.IndexStyle(::Type{<:PWM})
    return IndexLinear()
end

function Base.size(pwm::PWM)
    return size(pwm.data)
end

function Base.getindex(pwm::PWM, i::Integer)
    return getindex(pwm.data, i)
end

function Base.getindex(m::PWM{S}, s::S, j::Integer) where S<:Union{DNA,RNA}
    return m.data[index_nuc(s, j)]
end

function Base.show(io::IO, pwm::PWM{<:Union{DNA,RNA}})
    show_nuc_matrix(io, pwm)
end

function Base.show(io::IO, ::MIME"text/plain", m::PWM{<:Union{DNA,RNA}})
    show_nuc_matrix(io, m)
end

"""
    maxscore(pwm::PWM)

Return the maximum achievable score of `pwm`.
"""
function maxscore(pwm::PWM)
    if size(pwm, 2) == 0
        return zero(eltype(pwm))
    else
        return pwm.maxscore[1]
    end
end

# Make an accumulated maximum score vector.
function make_maxscore(pwm)
    len = size(pwm, 2)
    maxscore = Vector{eltype(pwm)}(undef, len)
    for j in len:-1:1
        s = pwm[1,j]
        for i in 2:size(pwm, 1)
            s = max(s, pwm[i,j])
        end
        if j == len
            maxscore[j] = s
        else
            maxscore[j] = s + maxscore[j+1]
        end
    end
    return maxscore
end

"""
    scoreat(seq::BioSequence, pwm::PWM, start::Integer)

Calculate the PWM score starting from `seq[start]`.
"""
function scoreat(seq::BioSequence, pwm::PWM, start::Integer)
    check_pwm(seq, pwm)
    pwmlen = size(pwm, 2)
    checkbounds(seq, start:start+pwmlen-1)
    score = zero(eltype(pwm))
    @inbounds for j in 0:pwmlen-1
        x = seq[start+j]
        score += iscertain(x) ? pwm[j<<2+trailing_zeros(x)+1] : zero(score)
    end
    return score
end

function Base.findfirst(pwm::PWM, seq::BioSequence, threshold::Real, start = 1, stop=lastindex(seq))
    if eltype(seq) == DNA || eltype(seq) == RNA
        return search_nuc(seq, start:stop, pwm, convert(eltype(pwm), threshold))
    else
        throw(ArgumentError("no search algorithm for '$(typeof(seq))'"))
    end
end

function check_pwm(seq, pwm::PWM{S}) where S <: Union{DNA,RNA}
    if size(pwm, 1) != 4
        throw(ArgumentError("matrix must have four rows"))
    elseif eltype(seq) != S
        throw(ArgumentError("sequence and PWM must have the same element type"))
    end
end

function search_nuc(seq::BioSequence, range::UnitRange{Int}, pwm::PWM{<:Union{DNA,RNA},S}, threshold::S) where S<:Real
    check_pwm(seq, pwm)
    checkbounds(seq, range)
    pwmlen = size(pwm, 2)
    for p in range.start:range.stop-pwmlen+1
        score = zero(eltype(pwm))
        @inbounds for j in 0:pwmlen-1
            if score + pwm.maxscore[j+1] < threshold
                break
            end
            x = seq[p+j]
            score += iscertain(x) ? pwm[j<<2+trailing_zeros(x)+1] : zero(score)
        end
        if score ≥ threshold
            return p
        end
    end
    return nothing
end

function index_nuc(s::Union{DNA,RNA}, j::Integer)
    if !iscertain(s)
        throw(ArgumentError(string("index symbol must be A, C, G or ", typeof(s) == DNA ? "T" : "U")))
    end
    return (j-1) << 2 + (trailing_zeros(s)) + 1
end

function show_nuc_matrix(io::IO, m::Union{PFM{S},PWM{S}}) where S<:Union{DNA,RNA}
    compact(x) = string(x ≥ 0 ? " " : "", sprint(show, x, context=:compact=>true))
    cells = hcat(['A', 'C', 'G', S == DNA ? 'T' : 'U'], compact.(m.data))
    width = maximum(length.(cells), dims=1)
    print(io, summary(m), ':')
    for i in 1:size(cells, 1)
        print(io, "\n ", rpad(cells[i,1], width[1]+1))
        for j in 2:size(cells, 2)
            print(io, rpad(cells[i,j], width[j]+1))
        end
    end
end
