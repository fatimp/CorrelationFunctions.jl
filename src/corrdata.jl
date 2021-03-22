const known_directions = [:x, :y, :z]

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

function mean(data :: CorrelationData, directions::Vector{Symbol} = known_directions)
    directions = check_directions(directions)
    return mapreduce(x -> data.success[x], +, directions) ./
        mapreduce(x -> data.total[x], +, directions)
end
