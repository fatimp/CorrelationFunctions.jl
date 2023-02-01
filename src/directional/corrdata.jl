"""
    CorrelationData

Structure returned by correlation functions (e.g. `l2`, `s2` or `c2`).

This structure holds correlation data computated along specified
directions. To access those data use the `mean` function.
"""
struct CorrelationData <: AbstractDict{AbstractDirection, Vector{Float64}}
    directions :: Vector{AbstractDirection}
    success    :: Dict{AbstractDirection, Vector{Float64}}
    total      :: Dict{AbstractDirection, Vector{Float64}}
end

function CorrelationData(len        :: Integer,
                         directions :: Vector{AbstractDirection})
    success = Dict(map(x -> x => zeros(len), directions))
    total   = Dict(map(x -> x => zeros(len), directions))
    return CorrelationData(directions, success, total)
end

# Copier
function CorrelationData(data :: CorrelationData)
    directions = data.directions
    success    = data.success
    total      = data.total

    return CorrelationData(directions,
                           Dict(x => success[x] for x in directions),
                           Dict(x => total[x]   for x in directions))
end

"""
    correlation_length(cd :: CorrelationData)

Return length of correlation vectors stored in `cd`.
"""
correlation_length(cd :: CorrelationData) = cd[cd |> keys |> first] |> length

# Undocumented AbstractDict interface
Base.keys(cd :: CorrelationData) = cd.directions
Base.haskey(cd :: CorrelationData, key) = key ∈ cd.directions
Base.length(cd :: CorrelationData) = length(cd.directions)

function Base.get(cd :: CorrelationData, key, default)
    if key ∈ cd.directions
        success = cd.success[key]
        total   = cd.total[key]
        return success ./ max.(total, 1)
    else
        return default
    end
end

function nextstate(cd :: CorrelationData, state)
    if state != nothing
        direction, next = state
        return direction => cd[direction], next
    else
        return nothing
    end
end

Base.iterate(cd :: CorrelationData) = nextstate(cd, iterate(cd.directions))
Base.iterate(cd :: CorrelationData, state) = nextstate(cd, iterate(cd.directions, state))

import StatsBase: mean
"""
    mean(data :: CorrelationData, directions :: Vector{AbstractDirection})
    mean(data :: CorrelationData)

Return mean of correlation function stored in `data` calculated along
directions stored in array `directions`. Elements of `directions` must
be members of `known_directions`. If `directions` contain only one
direction, values of correlation function computed along that
direction are returned. If `mean` called without `directions`
argument, mean of all computed directions is returned.
"""
function mean end

function mean(data :: CorrelationData, directions :: Vector{AbstractDirection})
    if !all(x -> x ∈ data.directions, directions)
        error("Correlation function is not computed for directions $directions")
    end

    success = mapreduce(x -> data.success[x],                  +, directions)
    # Once again, define correlation function f(a, x) = 0
    # for x > dimension of a.
    total   = mapreduce(x -> (y -> max(y, 1)).(data.total[x]), +, directions)

    return success ./ total
end

function mean(data :: CorrelationData)
    return mean(data, data.directions)
end
