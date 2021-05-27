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

    result = similar(img, Float64, result_size)

    CFMap(result, cf_type, img_size, zero_index)
end


function dir_from_map(cfmap::CFMap{T,N}, dir::Tuple{Vararg{Int,N}}) where T where N
    n = maximum(abs, dir)
    
    way = collect(0:n) / (n == 0 ? 1 : n)
    
    if cfmap.cf_type == :full
        ixs = map((z, d) -> round.(Int, way .* d .+ z), cfmap.zero_index, dir)

    elseif cfmap.cf_type == :central_symmetry
        if dir[end] < 0
            dir = .- dir
        end
        ixs = map((z, d) -> round.(Int, way .* d .+ z), cfmap.zero_index, dir)

    elseif cfmap.cf_type == :periodic_point_point
        ixs = map((z, d) -> round.(Int, way .* d .+ z), cfmap.zero_index, dir)
        for (ix, s) in zip(ixs, cfmap.img_size)
            for i in 1:length(ix)
                ix[i] += ix[i] < 1 ? s : 0
            end
        end
    else
        error()
    end
    
    cfmap.result[CartesianIndex.(ixs...)]
end


function restore_full_map(cfmap::CFMap{T,N}) where T where N
    res = cfmap.result

    if cfmap.cf_type == :full
        res

    elseif cfmap.cf_type == :central_symmetry
        n = size(res, N)
        reversed = reverse(res)
        [selectdim(reversed, N, 1:n - 1) res]

    elseif cfmap.cf_type == :periodic_point_point
        circshift(res, div.(size(res) .- 1, 2))
        
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
    elseif dir == :xy_main
        q = minimum(img_size) - 1
        (q, q)
    elseif dir == :xy_anti
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
        error(":$dir is not supported for 2D images")
    end
end


function dir_from_map(cfmap::CFMap, dir::Symbol)
    tdir = to_tuple_dir(cfmap.img_size, dir)
    dir_from_map(cfmap, tdir)
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
