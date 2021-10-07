# Convert neighborhood of a pixel into configuration index
local2conf(array :: AbstractArray{Bool, N}) where N =
    mapreduce(((idx, x),) -> x << (idx - 1), |, enumerate(array)) / (2^(2^N) - 1)

# TODO: write periodic version
function configuration_map(array :: AbstractArray{Bool})
    indices  = CartesianIndices(array)
    lidx = last(indices)
    uidx = oneunit(lidx)

    map(indices) do idx
        start = idx
        stop  = min(idx + uidx, lidx)
        slice = array[start:stop]
        local2conf(slice)
    end
end

function lc_crosscorr(array      :: AbstractArray{Bool};
                      directions :: Vector{Symbol} = array |> default_directions,
                      len        :: Integer = (array |> size |> minimum) ÷ 2,
                      periodic   :: Bool = false)
    return s2(configuration_map(array),
              SeparableIndicator(identity);
              len        = len,
              directions = directions,
              periodic   = periodic)
end
