"""
Structure returned by correlation functions (`l2`, `s2` and `c2`).

This structure holds correlation data computated along specified
directions. To access those data use the `mean` function.
"""
struct CorrelationData{T}
    directions :: Vector{Symbol}
    success    :: Dict{Symbol, Vector{T}}
    total      :: Dict{Symbol, Vector{T}}
end

function CorrelationData{T}(len          :: Integer,
                            directions   :: Vector{Symbol},
                            ndims        :: Integer) where T
    directions = check_directions(directions, ndims)
    success = Dict(map(x -> x => zeros(T, len), directions))
    total   = Dict(map(x -> x => zeros(T, len), directions))
    return CorrelationData{T}(directions, success, total)
end

function Base.getindex(x :: CorrelationData, i)
    directions = x.directions
    if i ∉ directions
        error("Requested direction $i is not among calculated directions $directions")
    end

    success = x.success[i]
    # Define correlation function f(a, x) = 0 for x > dimension of a
    total = (x -> max(x, 1)).(x.total[i])
    return success ./ total
end

function Base.show(io :: IO, x :: CorrelationData)
    directions = x.directions
    corr = reduce(hcat, x[direction] for direction in directions)
    pretty_table(io, corr, directions)
end

function Base.merge(d1 :: CorrelationData, d2 :: CorrelationData)
    return CorrelationData(unique(vcat(d1.directions, d2.directions)),
                           merge(d1.success, d2.success),
                           merge(d1.total,   d2.total))
end

import StatsBase: mean
"""
    mean(data :: CorrelationData, directions::Vector{Symbol})
    mean(data :: CorrelationData)

Return mean of correlation function stored in `data` calculated along
directions stored in array `directions`. Elements of `directions` must
be members of `known_directions`. If `directions` contain only one
direction, values of correlation function computed along that
direction are returned. If `mean` called without `directions`
argument, mean of all computed directions is returned.
"""
function mean end

function mean(data :: CorrelationData, directions::Vector{Symbol})
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

promote_accumulator_type(T :: Type) = (T <: Integer) ? Int : Float64
