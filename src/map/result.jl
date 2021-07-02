struct CFMap{T,N}
    result::T
    cf_type::Symbol
    img_size::Tuple{Vararg{Int,N}}
    zero_index::Tuple{Vararg{Int,N}}
end


function CFMap(img, cf_type)
    img_size = size(img)

    if cf_type == :full
        zero_index = img_size
        result_size = @. 2 * img_size - 1

    elseif cf_type == :central_symmetry
        zero_index = (img_size[1:end - 1]..., 1)
        result_size = (2 .* img_size[1:end - 1] .- 1 ..., img_size[end])

    elseif cf_type == :periodic_point_point
        zero_index = map(s -> 1, img_size)
        result_size = img_size

    else
        error("$cf_type is not supported")
    end

    result = similar(img, Float32, result_size)

    CFMap(result, cf_type, img_size, zero_index)
end


function dir_ixs(cfmap, dir)
    n = maximum(abs, dir)

    way = Float32.(collect(0:n) / (n == 0 ? 1 : n))

    if cfmap.cf_type == :central_symmetry && dir[end] < 0
        dir = .- dir
    end

    ixs = map((z, d) -> way .* d, cfmap.zero_index, dir)

    if cfmap.cf_type == :periodic_point_point
        for (ix, s) in zip(ixs, cfmap.img_size)
            out_ix = ix .< 0
            ix[out_ix] .+= s
        end
    end
    ixs
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
function dir_from_map(cfmap, dir)
    way = dir_ixs(cfmap, dir)
    itp = itp_map(cfmap)
    itp.(way...)
end


function surface_dirs(cfmap, dim, ispos)
    ixs = map(s -> -(s - 1):(s - 1), cfmap.img_size)
    map(ixs, 1:length(ixs)) do ix, d
        if d == dim
            if ispos
                ix[end]
            else
                ix[1]
            end
        else
            if d < dim
                ix[:]
            else
                ix[2:end-1]
            end
        end
    end
end


function all_dirs(cfmap)
    dirs = Vector{Tuple{Vararg{Int, ndims(cfmap.result)}}}()
    for i in 1:ndims(cfmap.result)
        for ispos in [true, false]
            
            s_dirs = surface_dirs(cfmap, i, ispos)
            for dir in Iterators.product(s_dirs...)
                push!(dirs, dir)
            end
        end
    end
    dirs
end


"""
    mean_dir(cfmap::CFMap)

Return averaged correlation function for r = 0, 1, 2, ...
"""
function mean_dir(cfmap::CFMap)
    data = cfmap.result
    
    m = cfmap.img_size .-1 |> norm |> ceil |> Int
    m += 1
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


"""
    restore_full_map(cfmap::CFMap)

Return centered full correlation function map.
"""
function restore_full_map(cfmap::CFMap{T,N}) where T where N
    res = cfmap.result

    if cfmap.cf_type == :full
        res

    elseif cfmap.cf_type == :central_symmetry
        n = size(res, N)
        reversed = reverse(res)
        [selectdim(reversed, N, 1:n - 1) res]

    elseif cfmap.cf_type == :periodic_point_point
        shiftmap(cfmap)[1]

    else
        error()
    end
end


function to_tuple_dir(img_size::Tuple{Int}, dir)
    if dir == :x
        (img_size[1] - 1,)
    else
        error(":$dir is not supported for 1D images")
    end
end


function to_tuple_dir(img_size::Tuple{Int,Int}, dir)
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


function to_tuple_dir(img_size::Tuple{Int,Int,Int}, dir)
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
    if cfmap.cf_type == :central_symmetry && mask[end]
        error("negative direction in last dimension is not supported for $(cfmap.cf_type)")
    elseif cfmap.cf_type == :periodic_point_point && any(mask)
        error("negative direction is not supported for $(cfmap.cf_type)")
    end

    ixs = map(cfmap.zero_index, cfmap.img_size, mask) do z, s, m
        if m
            z:-1:1
        else
            z:z + s - 1
        end
    end

    view(cfmap.result, ixs...)
end
