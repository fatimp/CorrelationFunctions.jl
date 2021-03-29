"""
Structure returned by correlation functions (`l2`, `s2` and `c2`).

This structure holds correlation data computated along specified
directions. To access those data use the `mean` function.
"""
struct CorrelationData
    directions :: Vector{Symbol}
    success    :: Dict{Symbol, Vector{Int}}
    total      :: Dict{Symbol, Vector{Int}}
end

function CorrelationData(len        :: Integer,
                         directions :: Vector{Symbol},
                         ndims      :: Integer)
    directions = check_directions(directions, ndims)
    success = Dict(map(x -> x => zeros(Int, len), directions))
    total   = Dict(map(x -> x => zeros(Int, len), directions))
    return CorrelationData(directions, success, total)
end

function Base.show(io :: IO, x :: CorrelationData)
    corr = mapreduce(hcat, x.directions) do direction
        # Define correlation function f(a, x) = 0 for x > dimension of a
        total = (x -> max(x, 1)).(x.total[direction])
        x.success[direction] ./ total
    end
    pretty_table(io, corr, x.directions)
end

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
