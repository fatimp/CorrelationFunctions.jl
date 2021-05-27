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
    x_L2_CUDA!(result, img, depth)

Compute L2 over first dim of `img` and *add* it to `result`
using CUDA.
"""
function x_L2_CUDA!(result, img, depth)
    indxs = CartesianIndices(size_slice(img, 1))
    nindxs = length(indxs)
    
    cu_index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    cu_stride = blockDim().x * gridDim().x

    for i in cu_index:cu_stride:nindxs
        indx = indxs[i]
        vres = view(result, :, indx)
        vimg = view(img, :, indx)
        base_L2!(vres, vimg, depth)
    end
    return
end


"""
when applied to CuArrays, performs computations on GPU
"""
function x_L2!(result::CuArray, img::CuArray, depth)
    result .= 0

    NTHREADS = 256
    N = reduce(*, size_slice(img, 1); init=1)
    
    NBLOCKS = div(N - 1, NTHREADS) + 1
    @cuda threads = NTHREADS blocks = NBLOCKS x_L2_CUDA!(result, img, depth)
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


function L2_positive_sides!(
    side_results,
    img::AbstractArray{<:Integer,N},
    side_depths,
    side_align_results,
    side_align_imgs;
    periodic
) where N
    for side in 1:N
        result = side_results[side]
        depth = side_depths[side]
        aligned_result = side_align_results[side]
        aligned_img = side_align_imgs[side]
        L2_side!(result, img, depth, side, aligned_result, aligned_img; periodic)
    end
end


function map_ix!(original_ixs, ray_ixs, img)
    N = length(original_ixs)
    for i in 1:N
        original_ixs[i] = []
        ray_ixs[i] = []
    end
    
    box = size(img)
    
    for orig_ix in CartesianIndices(img)
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
end


function restore!(result, side_results, original_ixs, ray_ixs)
    for (side_result, orig_ix, ray_ix) in zip(side_results, original_ixs, ray_ixs)
        result[orig_ix] .= side_result[ray_ix]
    end
end


struct Params_L2{Total,Result,AlignImg,N}
    # boundary conditions
    periodic::Bool
    # normalization
    total::Total

    # algorithm-specific

    side_results::Vector{Result}
    side_depths::Tuple{Vararg{Int,N}}
    side_align_imgs::Vector{AlignImg}
    side_align_results::Vector{Result}
    original_ixs::Vector{Vector{CartesianIndex{N}}}
    ray_ixs::Vector{Vector{CartesianIndex{N}}}
end


function l2(img::AbstractArray{<:Integer,N};
            periodic::Bool=true,
            depth::Int=img |> size |> maximum) where N
    
    side_sizes = [(size(img, d), size_slice(img, d)...) for d in 1:N]
    side_results = map(s -> similar(img, Int64, s), side_sizes)
    side_depths = map(s -> min(s, depth), size(img))

    if periodic
        side_align_img_sizes = map(s -> (2s[1], s[2:end]...), side_sizes)

        side_align_result_sizes = side_sizes
        
    else
        side_align_img_sizes = map(s -> (s[1], 2 .* s[2:end]...), side_sizes)

        side_align_result_sizes = side_align_img_sizes
    end
    side_align_imgs = map(s -> similar(img, Int8, s), side_align_img_sizes)
    side_align_results = map(s -> similar(img, Int, s...), side_align_result_sizes)
    
    
    original_ixs = Vector{Vector{CartesianIndex{N}}}(undef, N)
    ray_ixs = Vector{Vector{CartesianIndex{N}}}(undef, N)
    map_ix!(original_ixs, ray_ixs, img)
    
    
    p = Params_L2(
        periodic,
        cnt_total(img, periodic),
        side_results,
        side_depths,
        side_align_imgs,
        side_align_results,
        original_ixs,
        ray_ixs
    )
    cf_type = :central_symmetry
    p, cf_type
end


"""
    correllation_function!(res, img, params::Params_L2)

Compute L₂ correlation function in positive directions
"""
function correllation_function!(res, img, params::Params_L2)
    L2_positive_sides!(
        params.side_results,
        img,
        params.side_depths,
        params.side_align_results,
        params.side_align_imgs;
        params.periodic
    )
    restore!(res, params.side_results, params.original_ixs, params.ray_ixs)
    res ./= params.total
end
