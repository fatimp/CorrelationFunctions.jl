"""
~~~~
chord_length(array, phase;
             [directions = default_directions,]
             [nbins = 10,]
             [radius = 0.25,]
             [threshold = 0.3])
~~~~

Calculate chord length correlation function for a multiphase system
`array` and a phase `phase`.

Cord length function equals to probability of finding a chord whose
length is in the range `[l, l+δl]` and which lies entirely in the
phase `phase`. A chord is a line segment which touches boundary of the
phase `phase` with its ends.

This implementation bins chord lengths into `nbins` bins and returns
normalized histogram on collected data and the mean chord length in a
tuple.

`radius` and `threshold` parameters are used in edge detection
step. Do not change them if you don't know what you are doing.
"""
function chord_length(array      :: AbstractArray,
                      phase;
                      directions :: Vector{Symbol} = array |> ndims |> default_directions,
                      nbins      :: Integer       = 10,
                      radius     :: AbstractFloat = 0.25,
                      threshold  :: AbstractFloat = 0.3)
    # Select needed phase by applying indicator function to array
    ph = map(x -> x == phase, array)

    # Find a new "bogus" space for the edge
    edge_phase = maximum(array) + 1

    # List of chord lengths
    lengths = SLinkedList{Int}()

    for direction in directions
        # Edge detection along direction (like directional derivative)
        param = radius .* unit_vector(direction, ndims(array)) # Radii for Gauss filter
        blur = imfilter(ph, Kernel.gaussian(param))
        edge = abs.(blur - ph)

        # Apply threshold to the edge
        q = quantile(filter(x -> x != 0, edge), threshold)
        edge = map(x -> x > q ? edge_phase : 0, edge)

        # Combine the edge with the original picture
        combo = min.(edge + array, edge_phase)

        slicer = slice_generators(combo, Val(direction))
        for slice in slicer
            len = 0
            startonedge = false
            for x in slice
                # Edge found
                if x == edge_phase
                    len > 1 && startonedge && pushfirst!(lengths, len - 1)
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

    # List is not an AbstractVector, what a shame!
    len_array = [x for x in lengths]
    hist = fit(Histogram, len_array; nbins = nbins)
    mean_len = sum(len_array) ./ length(len_array)
    return normalize(hist; mode = :probability), mean_len
end
