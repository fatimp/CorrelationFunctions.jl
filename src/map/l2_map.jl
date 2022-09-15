"""
    base_L2!(result, img, depth=length(result))

Compute L2 over vector `img` and *add* it to `result`.

Simple and fast L2 algorithm,
`O(n), n = length(img)`.

Meaningfull values for elements of `img`:
* -1 -- not a part of the image
* 0 -- void
* 1 -- solid, inside original boundaries
* 2 -- solid, outside original boundaries
"""
function base_L2!(result, img, depth=length(result))
    len = 0
    len_in_original = 0

    for (i, x) in enumerate(img)
        if x >= 1
            len += 1
            if x == 1
                len_in_original += 1
            end
        end

        if x < 1 || i == length(img)
            for r in 1:min(depth, len)
                result[r] += min(len - r + 1, len_in_original)
            end
            len = 0
            len_in_original = 0
        end
    end
end

function align!(aimg, img, side, ray, periodic)
    n = size(img, side)
    for i in 1:n
        aslice = selectdim(aimg, 1, i)
        slice = selectdim(img, side, i)

        j = round.(Int, (i - 1) .* ray ./ (n - 1))

        if periodic
            circshift!(aslice, slice, .-j)

            aslice2 = selectdim(aimg, 1, i + n)
            k = round.(Int, (n + i - 1) .* ray ./ (n - 1))
            circshift!(aslice2, slice, .-k)
            aslice2 .*= 2
        else
            ix = CartesianIndex(size(slice)) .+ CartesianIndices(slice) .- CartesianIndex(j)
            aslice .= -1
            aslice[ix] .= slice
        end
    end

    return aimg
end

function L2_side(img, side, periodic)
    side_size = (size(img, side), size_slice(img, side)...)
    side_align_img_size = periodic ?
        (2side_size[1], side_size[2:end]...) :
        (side_size[1], 2 .* side_size[2:end]...)
    side_align_result_size = periodic ? side_size : side_align_img_size

    side_result = zeros(Int, side_size)
    aligned_img = zeros(Int, side_align_img_size)

    ray_ixs = CartesianIndices(size_slice(img, side))

    for ray_ix in ray_ixs
        ray_result = view(side_result, :, ray_ix)
        ray_projection = Tuple(ray_ix) .- 1

        aligned_img = align!(aligned_img, img, side, ray_projection, periodic)
        indices = CartesianIndices(size_slice(aligned_img, 1))
        aligned_result = zeros(Int, side_align_result_size)
        for index in indices
            data_x = view(aligned_img, :, index)
            result_x = view(aligned_result, :, index)
            base_L2!(result_x, data_x, size(img, side))
        end
        sum!(ray_result, aligned_result)
    end

    return side_result
end


function L2_positive_sides(img :: AbstractArray, original_ixs, ray_ixs,
                           periodic)
    result = zeros(Float32, size(img))

    for side in 1:ndims(img)
        orig_ix = original_ixs[side]
        ray_ix = ray_ixs[side]
        side_result = L2_side(img, side, periodic)
        result[orig_ix] .= side_result[ray_ix]
    end

    return result
end

function map_ix(indices :: CartesianIndices{N}) where N
    box = size(indices)

    original_ixs = [CartesianIndex{N}[] for _ in 1:N]
    ray_ixs = [CartesianIndex{N}[] for _ in 1:N]

    for orig_ix in indices
        coords = Tuple(orig_ix) .- 1
        side = argmax(coords ./ box)

        n = box[side] - 1
        q = coords[side]

        if q == 0
            rayend = box
        else
            rayend = round.(Int, coords ./ q .* n) .+ 1
        end

        mask = [1:side - 1; side + 1:N]
        rayend = rayend[mask]

        ray_ix = CartesianIndex(q + 1, rayend...)
        push!(original_ixs[side], orig_ix)
        push!(ray_ixs[side], ray_ix)
    end

    return original_ixs, ray_ixs
end

"""
    correllation_function!(res, img, params::Params_L2)

Compute L₂ correlation function in positive directions
"""
function correllation_function(img, orig_ixs, ray_ixs, periodic)
    res = L2_positive_sides(img, orig_ixs, ray_ixs, periodic)
    return normalize_result(res, periodic)
end

bool_mask(x, n) = digits(Bool, x; base = 2, pad = n)
bool_iter(n) = (bool_mask(x, n) for x in 0:(2^n-1))

@doc raw"""
    l2(image, phase; periodic = false)

Calculate $L_2$ (lineal path) correlation function for the phase
`phase` in an N-dimensional image (up to 3 dimensions) and return a
`CFMap` object.

# Examples
```jldoctest
julia> l2([1 0; 0 1], 1; periodic=true).result
3×2 Matrix{Float64}:
 0.0  0.5
 0.5  0.0
 0.0  0.5
```
"""
function l2(image :: AbstractArray, phase; periodic=false)
    phase_array = image .== phase
    result = CFMap(phase_array)
    orig_ixs, ray_ixs = image |> CartesianIndices |> map_ix

    for mask in bool_iter(ndims(phase_array))
        if mask[end]
            continue
        end

        mirror_img = mirror(phase_array, mask)
        mirror_result = correllation_function(mirror_img, orig_ixs, ray_ixs, periodic)

        v_result = cut_cfmap(result, mask)
        v_result .= mirror_result
    end

    return result
end
