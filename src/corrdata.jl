"""
`known_directions` array holds directions along which correlation
functions can be computed. Currently computations can be done in the
following way:

* Along three axes whose directions are designated as `:x`, `:y`
  and `:z`.
* Along diagonals of 2D slices of the input array. These 2D slices
  can have the equations `x = C`, `y = C` and `z = C` where C is some
  constant and the corresponding directions are designated as
  `:yz_main`, `:xz_main` and `:xy_main`. Diagonals used in this case
  are parallel to the main diagonals of 2D slices.
* Along diagonals of the same 2D slices, but in this case the
  diagonals are parallel to the antidiagonal of the slices. These
  directions are designated as `:yz_anti`, `:xz_anti` or `:xy_anti`.
"""
const known_directions = [:x,       :y,       :z,
                          :yz_main, :xz_main, :xy_main,
                          :yz_anti, :xz_anti, :xy_anti]

"""
`default_directions` is a subset of `known_directions` which is used
as a default value of `directions` argument for `l2`, `s2` and `c2`
functions.
"""
const default_directions = [:x, :y, :z]

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
    directions = check_directions(directions)
    return mapreduce(x -> data.success[x], +, directions) ./
        mapreduce(x -> data.total[x], +, directions)
end

function mean(data :: CorrelationData)
    return mean(data, data.directions)
end
