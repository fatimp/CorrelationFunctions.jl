struct CFMap{T,N}
    result::T
    cf_type::Symbol
    img_size::Tuple{Vararg{Int,N}}
    zero_index::Tuple{Vararg{Int,N}}
end


function CFMap(img, cf_type)
    @assert cf_type == :central_symmetry
    img_size = size(img)
    zero_index = (img_size[1:end - 1]..., 1)
    result_size = (2 .* img_size[1:end - 1] .- 1 ..., img_size[end])
    result = similar(img, Float32, result_size)

    return CFMap(result, cf_type, img_size, zero_index)
end


function dir_ixs(cfmap, dir)
    @assert cfmap.cf_type == :central_symmetry

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


function shiftmap(cfmap::CFMap)
    res = cfmap.result
    ns = size(res)
    ms = ns .÷ 2
    raw_ixs = map(m -> -m:m, ms)
    ixs = map((rix, n) -> mod.(rix, n) .+ 1, raw_ixs, ns)
    res[ixs...], raw_ixs
end


function itp_map(cfmap)
    if cfmap.cf_type == :periodic_point_point
        shiftres, xs = shiftmap(cfmap)
        itp = interpolate(xs, shiftres, Gridded(Linear()))
        extrapolate(itp, Periodic())
    else
        xs = map((ix, z) -> ix .- z, axes(cfmap.result), cfmap.zero_index)
        itp = interpolate(xs, cfmap.result, Gridded(Linear()))
    end
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

"""
    mean_dir(cfmap::CFMap)

Return averaged correlation function for r = 0, 1, 2, ...
"""
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
    @assert cfmap.cf_type == :central_symmetry

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
    @assert cfmap.cf_type == :central_symmetry
    if mask[end]
        error("negative direction in last dimension is not supported for $(cfmap.cf_type)")
    end

    ixs = map(cfmap.zero_index, cfmap.img_size, mask) do z, s, m
        m ? (z:-1:1) : (z:z + s - 1)
    end

    return view(cfmap.result, ixs...)
end
