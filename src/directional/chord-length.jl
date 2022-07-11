"""
Structure returned by `chord_length` function. To access chord length
histogram use `hist` field.
"""
struct ChordLengthInfo
    μ :: Float64
    σ :: Float64
    hist
end

function Base.show(io :: IO, ::MIME"text/plain", x :: ChordLengthInfo)
    print(io, "Chord length info (mean = $(x.μ), std = $(x.σ))")
end

"""
    chord_length(array, phase[; directions,] nbins = 10)

Calculate chord length correlation function for one-, two- or
three-dimensional multiphase system.

Cord length function `p(x)` equals to probability of finding a chord
whose length is in the range `[x, x+δx]` and which lies entirely in the
phase `phase`. A chord is a line segment which touches the boundary of
a same-phase cluster with its ends.

This implementation bins chord lengths into `nbins` bins and returns
normalized histogram on collected data along with mean chord length
and standard deviation.

# Examples
```jldoctest
julia> chord_length([1, 0, 0, 0, 0, 1, 0, 1], 0)
Chord length info (mean = 2.5, std = 2.1213203435596424)
```

For a list of possible dimensions, see also: [`direction1Dp`](@ref),
[`direction2Dp`](@ref), [`direction3Dp`](@ref).
"""
function chord_length(array      :: AbstractArray,
                      phase;
                      directions :: Vector{Symbol} = array |> default_directions,
                      nbins      :: Integer       = 10)
    # Select needed phase by applying indicator function to array
    ph = map(x -> x == phase, array)

    # Find a new "bogus" space for the edge
    edge_phase = maximum(array) + 1

    # Arary of chord lengths
    lengths = Vector{Int}(undef, 0)

    # Poor man's edge detection
    dist = ph |> feature_transform |> distance_transform
    edge = edge_phase * map(x -> x == 1.0, dist)

    # Combine the edge with the original picture
    combo = min.(edge + array, edge_phase)
    for direction in directions
        slicer = slice_generators(combo, false, Val(direction))
        for slice in slicer
            len = 0
            startonedge = false
            for x in slice
                # Edge found
                if x == edge_phase
                    len > 1 && startonedge && push!(lengths, len - 1)
                    len = 0
                    startonedge = true
                # Not our phase
                elseif x != phase
                    startonedge = false
                end
                len += 1
            end
        end
    end

    hist = fit(Histogram, lengths; nbins = nbins)
    return ChordLengthInfo(mean(lengths),
                           std(lengths),
                           normalize(hist; mode = :probability))
end
