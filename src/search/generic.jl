
function Base.findall(pat, seq::BioSequence, start::Int=1, stop::Int=length(seq))

    results = Vector{UnitRange{Int}}()

    while start <= stop

        r = findfirst(pat, seq, start, stop) # Note returns the range of the first occurrence.

        if r === nothing
            break
        end

        push!(results, r)

        start = last(r) + 1
    end

    return results
end

function Base.findall(pat, seq::BioSequence, range::UnitRange{<:Integer} = 1:length(seq))

    start = first(range)
    stop = last(range)

    return findall(pat, seq, start, stop)
end

function Base.findall(pat, seq::BioSequence, ranges::AbstractVector{R}) where {R<:UnitRange{<:Integer}}

    results = Vector{R}()

    for range in ranges
        results = vcat(results, findall(pat, seq, range))
    end

    return results
end
