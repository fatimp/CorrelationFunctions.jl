"""
known_directions array holds directions along which correlation
functions can be computated. Currently computations can be done along
three array axes and directions are designated as `:x`, `:y` and
`:z`.
"""
const known_directions = [:x, :y, :z]

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

function check_directions(directions :: Vector{Symbol})
    directions = unique(directions)
    if !all(x -> x ∈ known_directions, directions)
        error("Directions must be members of $known_directions")
    end

    return directions
end

function CorrelationData(len :: Integer, directions::Vector{Symbol})
    directions = check_directions(directions)
    success = Dict(map(x -> x => zeros(Int, len), directions))
    total   = Dict(map(x -> x => zeros(Int, len), directions))
    return CorrelationData(directions, success, total)
end

function Base.show(io :: IO, x :: CorrelationData)
    corr = mapreduce(direction -> x.success[direction] ./ x.total[direction], hcat, x.directions)
    pretty_table(io, corr, x.directions)
end

"""
    mean(data :: CorrelationData, directions::Vector{Symbol} = known_directions)

Return values of correlation function stored in `data` averaged from
directions specified in the array `directions`. Elements of
`directions` must be members of `known_directions`. If `directions`
contain only one direction, values of correlation function computed
along that direction are returned.
"""
function mean(data :: CorrelationData, directions::Vector{Symbol} = known_directions)
    directions = check_directions(directions)
    return mapreduce(x -> data.success[x], +, directions) ./
        mapreduce(x -> data.total[x], +, directions)
end
