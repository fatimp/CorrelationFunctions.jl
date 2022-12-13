"""
    dir_from_map(m, dir)

Extract a direction from a correlation map. `direction` can be either
`:x`, `:y`, `:z`, `:xy` or `:yx`.
"""
function dir_from_map(m::AbstractArray, direction; periodic=false)
    if direction == :x
        q = periodic ? size(m, 1) : size(m, 1) ÷ 2 + 1
        ixs = CartesianIndex.(1:q, 1, 1)
    elseif direction == :y
        q = periodic ? size(m, 2) : size(m, 2) ÷ 2 + 1
        ixs = CartesianIndex.(1, 1:q, 1)
    elseif direction == :z
        q = periodic ? size(m, 3) : size(m, 3) ÷ 2 + 1
        ixs = CartesianIndex.(1, 1, 1:q)
    elseif direction == :xy
        a = minimum(size(m)[1:2])
        q = periodic ? a : a ÷ 2 + 1
        ixs = CartesianIndex.(1:q, 1:q, 1)
    elseif direction == :yx
        a = minimum(size(m)[1:2])
        q = periodic ? a : a ÷ 2 + 1
        b = 1:q
        c = @. mod(1 - b, a) + 1
        ixs = CartesianIndex.(c, b, 1)
    end
    m[ixs]
end

"""
    average_directions(cfmap; len = (cfmap |> size |> minimum) ÷ 2)

Average correlation map `cfmap` over all directions. The result is a
vector of length `len` with indices being equal to correlation length + 1.
"""
function average_directions(cfmap :: AbstractArray{T};
                            len   :: Integer = (cfmap |> size |> minimum) ÷ 2) where T
    counter = zeros(Int, len)
    accum   = zeros(T,   len)

    for idx in CartesianIndices(cfmap)
        dist = round(Int, ((Tuple(idx) .- 1) .^ 2) |> sum |> sqrt) + 1
        if dist <= len
            accum[dist] += cfmap[idx]
            counter[dist] += 1
        end
    end

    return accum ./ counter
end

#====================#
# Fuck you Julia! Why maybe_upload_to_gpu() causes "type instability"
# here but not in extract_edges()?! The worst language ever!
function cnt_total(c :: AbstractArray)
    return map(axes(c), size(c)) do ix, es
        s = (es + 1) ÷ 2
        cnt = @. s - abs(ix - s)
        cnt |> ifftshift
    end
end

function cnt_total(c :: CuArray)
    return map(axes(c), size(c)) do ix, es
        s = (es + 1) ÷ 2
        cnt = @. s - abs(ix - s)
        cnt |> ifftshift |> CuArray
    end
end
#====================#

cnt_total_reshaped(c :: AbstractArray{T, 1}) where T = cnt_total(c)
cnt_total_reshaped(c :: AbstractArray{T, 2}) where T =
    let (t1, t2) = cnt_total(c);
        (reshape(t1, :, 1), reshape(t2, 1, :))
    end
cnt_total_reshaped(c :: AbstractArray{T, 3}) where T =
    let (t1, t2, t3) = cnt_total(c);
        (reshape(t1, :, 1, 1), reshape(t2, 1, :, 1), reshape(t3, 1, 1, :))
    end

function normalize_result(result   :: AbstractArray,
                          periodic :: Bool)
    local norm

    if periodic
        norm = result / length(result)
    else
        total = cnt_total_reshaped(result)
        norm = reduce(./, total; init = result)
    end

    return norm
end

function zeropad(image)
    s = size(image)
    s = (2 .* s) .- 1

    padded = similar(image, s)
    padded .= 0
    padded[axes(image)...] .= image
    return padded
end
