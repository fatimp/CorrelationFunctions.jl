"""
`default_directions_3d` is a subset of all possible direction which is
used as a default value of `directions` argument for `l2`, `s2` and `c2`
functions when using with 3D array.
"""
const default_directions_3d = [:x, :y, :z]

"""
`default_directions_3d` is a subset of all possible direction which is
used as a default value of `directions` argument for `l2`, `s2` and `c2`
functions when using with 2D array.
"""
const default_directions_2d = [:x, :y]

"""
    default_directions(ndims)

Get default direction in which correlation functions are calculated
for array of this number of dimensions.
"""
function default_directions end

default_directions(x :: Integer) = x |> Val |> default_directions
default_directions(:: Val{3}) = default_directions_3d
default_directions(:: Val{2}) = default_directions_2d

"""
    direction2Dp(sym)

Return true is `sym` is 2D direction or false otherwise. Known
directions are:

* `x` and `y`. Correlation functions computed in those directions are
   computed on slices taken along those axes.
* `xy_main`. Correlation functions are computed in the direction
   parallel to the main diagonal on the array.
* `xy_anti`. Correlation functions are computed in the direction
   parallel to the antidiagonal on the array.
"""
function direction2Dp end

direction2Dp(x :: Symbol) = x |> Val |> direction2Dp
direction2Dp(  :: Any) = false

macro define_2d_direction(sym)
    return :(direction2Dp(:: Val{$sym}) = true)
end

@define_2d_direction :x
@define_2d_direction :y
@define_2d_direction :xy_main
@define_2d_direction :xy_anti

"""
    direction3Dp(sym)

Return true is `sym` is 3D direction or false otherwise. Known
directions are:

* `x`, `y` and `z`. Correlation functions computed in those directions
   are computed on slices taken along those axes.
* `xy_main`, `xz_main` and `yz_main`. To compute correlation functions
   in those directions 2D planes with equations `z = const`, `y = const`
   and `x = const` respectively are cut from the input data and
   computations are done in direction parallel to the main diagonal of
   those slices.
* `xy_anti`, `xz_anti` and `yz_anti`. The same as above, only
   computations are done in directions parallel to the antidiagonal of
   slices.
"""
function direction3Dp end

direction3Dp(x :: Symbol) = x |> Val |> direction3Dp
direction3Dp(  :: Any) = false

macro define_3d_direction(sym)
    return :(direction3Dp(:: Val{$sym}) = true)
end

@define_3d_direction :x
@define_3d_direction :y
@define_3d_direction :z
@define_3d_direction :xy_main
@define_3d_direction :yz_main
@define_3d_direction :xz_main
@define_3d_direction :xy_anti
@define_3d_direction :yz_anti
@define_3d_direction :xz_anti

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

function check_directions(directions :: Vector{Symbol},
                          ndims      :: Integer)
    if ndims == 3
        predicate = direction3Dp
    elseif ndims == 2
        predicate = direction2Dp
    else
        error("Wrong number of dimensions")
    end

    directions = unique(directions)

    if !all(predicate, directions)
        error("Unknown dimensions found. See documentation to
        `direction3Dp` or `direction2Dp`")
    end

    return directions
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
