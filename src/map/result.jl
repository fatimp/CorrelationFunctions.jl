struct CFMap{T,N}
    result     :: T
    img_size   :: Tuple{Vararg{Int,N}}
    zero_index :: Tuple{Vararg{Int,N}}
end


function CFMap(img)
    img_size = size(img)
    zero_index = (img_size[1:end - 1]..., 1)
    result_size = (2 .* img_size[1:end - 1] .- 1 ..., img_size[end])
    result = similar(img, Float32, result_size)

    return CFMap(result, img_size, zero_index)
end


function dir_ixs(cfmap, dir)
    n = maximum(abs, dir)
    way = Float32.(collect(0:n) / (n == 0 ? 1 : n))

    if dir[end] < 0
        dir = .- dir
    end

    return map((z, d) -> way .* d, cfmap.zero_index, dir)
end


function dir_ixs(cfmap, dir::Symbol)
    tdir = to_tuple_dir(cfmap.img_size, dir)
    dir_ixs(cfmap, tdir)
end

function itp_map(cfmap)
    xs = map((ix, z) -> ix .- z, axes(cfmap.result), cfmap.zero_index)
    return interpolate(xs, cfmap.result, Gridded(Linear()))
end


"""
    dir_from_map(cfmap, dir)

Return CF on one direcion using interpolations.

`CF.(r...)`, where `r = Map.dir_ixs(cfmap, dir)`. 
"""
function dir_from_map(cfmap::CFMap, dir)
    way = dir_ixs(cfmap, dir)
    itp = itp_map(cfmap)
    itp.(way...)
end

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
    mean_dir(cfmap)

Average correlation map `cfmap` over all directions. The result is a
vector with indices being equal to correlation length + 1.
"""
function mean_dir end

function mean_dir(cfmap::CFMap)
    data = cfmap.result
    
    m = cfmap.img_size |> norm |> floor |> Int
    m += ndims(data) > 1
    
    cfdir = zeros(Float32, m)
    weights = zeros(Float32, m)
    
    for ix in CartesianIndices(data)
        tix = Tuple(ix)
        r = norm(tix .- cfmap.zero_index)
        
        q = floor(Int, r) + 1
        α = q - r

        cfdir[q] += data[ix] * α
        cfdir[q + 1] += data[ix] * (1 - α)
        weights[q] += α
        weights[q + 1] += (1 - α)
    end
    
    cfdir ./= weights
end

function mean_dir(cfmap :: AbstractArray{T}) where T
    num  = Dict{Int, Int}()
    vals = Dict{Int, Vector{T}}()
    for idx in CartesianIndices(cfmap)
        dist = ((Tuple(idx) .- 1) .^ 2) |> sum |> sqrt |> floor |> Int
        v = get(vals, dist, T[])
        push!(v, cfmap[idx])
        num[dist]  = get(num, dist, 0) + 1
        vals[dist] = v
    end

    len = num |> keys |> maximum
    avg = Vector{Float64}(undef, len)

    for idx in 1:len
        n = get(num, idx-1, 0)
        avg[idx] = iszero(n) ? 0 : (sum(vals[idx-1]) / n)
    end

    return avg
end

"""
    restore_full_map(cfmap::CFMap)

Return centered full correlation function map.
"""
function restore_full_map(cfmap::CFMap{T,N}) where {T, N}
    res = cfmap.result
    n = size(res, N)
    reversed = reverse(res)
    return [selectdim(reversed, N, 1:n - 1) res]
end


function to_tuple_dir(img_size :: NTuple{1, Int}, dir)
    if dir == :x
        (img_size[1] - 1,)
    else
        error(":$dir is not supported for 1D images")
    end
end


function to_tuple_dir(img_size :: NTuple{2, Int}, dir)
    if dir == :x
        (img_size[1] - 1, 0)
    elseif dir == :y
        (0, img_size[2] - 1)
    elseif dir == :xy
        q = minimum(img_size) - 1
        (q, q)
    elseif dir == :yx
        q = minimum(img_size) - 1
        (-q, q)
    else
        error(":$dir is not supported for 2D images")
    end
end


function to_tuple_dir(img_size :: NTuple{3, Int}, dir)
    if dir == :x
        (img_size[1] - 1, 0, 0)
    elseif dir == :y
        (0, img_size[2] - 1, 0)
    elseif dir == :z
        (0, 0, img_size[3] - 1)
    else
        error(":$dir is not supported for 3D images")
    end
end


function cut_cfmap(cfmap, mask)
    if mask[end]
        error("negative direction in last dimension is not supported")
    end

    ixs = map(cfmap.zero_index, cfmap.img_size, mask) do z, s, m
        m ? (z:-1:1) : (z:z + s - 1)
    end

    return view(cfmap.result, ixs...)
end
