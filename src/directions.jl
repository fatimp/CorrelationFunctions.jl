"""
    direction1Dp(sym)

Return true is `sym` is 1D direction or false otherwise. Only returns
`true` if `sym` is `:x`.
"""
function direction1Dp end

"""
    direction2Dp(sym)

Return true is `sym` is 2D direction or false otherwise. Known
directions are:

* `:x` and `:y`. Correlation functions are computed along unit vectors
   `(0, 1)` and `(1, 0)` respectively.
* `:xy_main` and `:xy_anti`. Correlation functions are computed in
   diagonal directions `(1, 1)` and `(-1, 1)` respectively.
"""
function direction2Dp end

"""
    direction3Dp(sym)

Return true is `sym` is 3D direction or false otherwise. Known
directions are:

* `:x`, `:y` and `:z`. Correlation functions are computed along unit
   vectors `(1, 0, 0)`, `(0, 1, 0)` and `(0, 0, 1)` respectively.
* `:xy_main`, `:xz_main` and `:yz_main`. Correlation functions are
   computed in diagonal directions `(1, 1, 0)`, `(1, 0, 1)` and
   `(0, 1, 1)` respectively.
* `:xy_anti`, `:xz_anti` and `:yz_anti`. Correlation functions are
   computed in diagonal directions `(-1, 1, 0)`, `(-1, 0, 1)` and
   `(0, -1, 1)` respectively.
* From `:diag1` to `:diag4`. Corresponding directions are `(1, 1, 1)`,
   `(-1, 1, 1)`, `(1, -1, 1)` and `(1, 1, -1)`.
"""
function direction3Dp end

macro def_direction_predicate(name)
    exp = 
        quote
            $name(x :: Symbol) = x |> Val |> $name
            $name(  :: Any) = false
        end
    return esc(exp)
end

macro def_direction(predicate, sym)
    return :($predicate(:: Val{$sym}) = true)
end

@def_direction_predicate direction1Dp
@def_direction_predicate direction2Dp
@def_direction_predicate direction3Dp

@def_direction(direction1Dp, :x)

@def_direction(direction2Dp, :x)
@def_direction(direction2Dp, :y)
@def_direction(direction2Dp, :xy_main)
@def_direction(direction2Dp, :xy_anti)

@def_direction(direction3Dp, :x)
@def_direction(direction3Dp, :y)
@def_direction(direction3Dp, :z)
@def_direction(direction3Dp, :xy_main)
@def_direction(direction3Dp, :yz_main)
@def_direction(direction3Dp, :xz_main)
@def_direction(direction3Dp, :xy_anti)
@def_direction(direction3Dp, :yz_anti)
@def_direction(direction3Dp, :xz_anti)

# True diagonals (TODO: better names?)
@def_direction(direction3Dp, :diag1)
@def_direction(direction3Dp, :diag2)
@def_direction(direction3Dp, :diag3)
@def_direction(direction3Dp, :diag4)

"""
    default_directions(ndims)

Get default direction in which correlation functions are calculated
for array of this number of dimensions.
"""
function default_directions end

default_directions(x :: Integer) = x |> Val |> default_directions
default_directions(:: Val{3}) = [:x, :y, :z]
default_directions(:: Val{2}) = [:x, :y]
default_directions(:: Val{1}) = [:x]

"""
    direction_predicate(ndims)

Get direction predicate for specified number of dimensions `ndims`.
"""
function direction_predicate end

direction_predicate(x :: Integer) = x |> Val |> direction_predicate
direction_predicate(:: Val{3}) = direction3Dp
direction_predicate(:: Val{2}) = direction2Dp
direction_predicate(:: Val{1}) = direction1Dp
direction_predicate(:: Any) = error("Wrong number of dimensions")

function check_directions(directions :: Vector{Symbol},
                          ndims      :: Integer)
    predicate = direction_predicate(ndims)
    directions = unique(directions)

    if !all(predicate, directions)
        error("Unknown dimensions found.")
    end

    return directions
end
