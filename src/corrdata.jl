const known_directions = [:x, :y, :z]

struct CorrelationData
    directions :: Vector{Symbol}
    success    :: Dict{Symbol, Vector{Int}}
    total      :: Dict{Symbol, Vector{Int}}
end

function check_directions(directions)
    dirs :: Vector{Symbol} = unique(directions)
    if !all(x -> x ∈ known_directions, dirs)
        error("Directions must be members of $known_directions")
    end

    return dirs
end

function CorrelationData(len :: Integer, directions::Symbol...)
    dirs = check_directions(directions)
    success = Dict(map(x -> x => zeros(Int, len), dirs))
    total   = Dict(map(x -> x => zeros(Int, len), dirs))
    return CorrelationData(dirs, success, total)
end

function Base.show(io :: IO, x :: CorrelationData)
    corr = mapreduce(direction -> x.success[direction] ./ x.total[direction], hcat, x.directions)
    pretty_table(io, corr, x.directions)
end

function mean(data :: CorrelationData, directions::Symbol...)
    dirs = check_directions(directions)
    return mapreduce(x -> data.success[x], +, dirs) ./ mapreduce(x -> data.total[x], +, dirs)
end

mean(data :: CorrelationData) = mean(data, data.directions...)
