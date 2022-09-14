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


"""
    x_L2!(result, img, depth)

Compute L2 over first dim of `img` and *add* it to `result`.

Simple and fast single-CPU algorithm: `O(N)`,
where `N = length(img)`.
"""
function x_L2!(result, img, depth)
    result .= 0

    indxs = CartesianIndices(size_slice(img, 1))
    for indx in indxs
        data_x = view(img, :, indx)
        result_x = view(result, :, indx)
        base_L2!(result_x, data_x, depth)
    end
    return
end

"""
    align!(aimg, img, side, ray; periodic=true)

"""
function align!(aimg, img, side, ray; periodic=true)
# TODO well documented align
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
    aimg
end

"""
    align!(aimg, img, side, ray; periodic=true)

Simple align for 1D
"""
function align!(aimg::AbstractVector, img::AbstractVector, side, ray; periodic=true)
    if periodic
        n = length(img)
        aimg[1:n] .= img
        aimg[n + 1:end] .= 2 .* img
    else
        aimg .= img
    end
    aimg
end


function L2_side!(
    side_result,
    img,
    depth,
    side,
    aligned_result,
    aligned_img;
    periodic::Bool
)
    ray_ixs = CartesianIndices(size_slice(img, side))

    for ray_ix in ray_ixs
        ray_result = view(side_result, :, ray_ix)
        ray_projection = Tuple(ray_ix) .- 1

        align!(aligned_img, img, side, ray_projection; periodic)
        x_L2!(aligned_result, aligned_img, depth)
        sum!(ray_result, aligned_result)
    end
end


function L2_positive_sides(
    img :: AbstractArray,
    original_ixs,
    ray_ixs;
    periodic
)
    result = zeros(Float32, size(img))

    for side in 1:ndims(img)
        side_size = (size(img, side), size_slice(img, side)...)


        side_align_img_size = periodic ? (2side_size[1], side_size[2:end]...) :
            (side_size[1], 2 .* side_size[2:end]...)
        side_align_result_size = periodic ? side_size : side_align_img_size

        side_result = zeros(Int, side_size)
        depth = size(img, side)
        orig_ix = original_ixs[side]
        ray_ix = ray_ixs[side]
        aligned_result = zeros(Int, side_align_result_size)
        aligned_img = zeros(Int, side_align_img_size)
        L2_side!(side_result, img, depth, side, aligned_result, aligned_img; periodic)
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

struct Params_L2{N}
    # boundary conditions
    periodic::Bool

    # algorithm-specific

    original_ixs::Vector{Vector{CartesianIndex{N}}}
    ray_ixs::Vector{Vector{CartesianIndex{N}}}
end


function Params_L2(img::AbstractArray{<:Integer,N};
                   periodic::Bool=true) where N
    original_ixs = Vector{Vector{CartesianIndex{N}}}(undef, N)
    ray_ixs = Vector{Vector{CartesianIndex{N}}}(undef, N)
    original_ixs, ray_ixs = img |> CartesianIndices |> map_ix

    return Params_L2(periodic,
                     original_ixs,
                     ray_ixs)
end


"""
    correllation_function!(res, img, params::Params_L2)

Compute L₂ correlation function in positive directions
"""
function correllation_function(img, params::Params_L2)
    res = L2_positive_sides(
        img,
        params.original_ixs, params.ray_ixs;
        params.periodic
    )
    return normalize_result(res, params.periodic)
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
    cf_params = Params_L2(phase_array; periodic = periodic)
    result = CFMap(phase_array)

    for mask in bool_iter(ndims(phase_array))
        if mask[end]
            continue
        end

        mirror_img = mirror(phase_array, mask)
        mirror_result = correllation_function(mirror_img, cf_params)

        v_result = cut_cfmap(result, mask)
        v_result .= mirror_result
    end

    return result
end
