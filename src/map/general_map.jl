function mirror(image, mask)
    ix = map((i, m) -> m ? reverse(i) : i, axes(image), mask)
    
    view(image, ix...)
end


function index_res(ix, shift, mask, rmask)
    if !rmask
        ix
    elseif mask
        reverse(ix)
    else
        return ix .+ shift
    end
end


function maskres(res, image, mask, rmask)
    # TODO резать на основе размерностей res, а не image
    ix = map(index_res, axes(image), size(image), mask, rmask)
    view(res, ix...)
end


"""
    Params_map

Store parameters and inner pre-allocated matrices.
"""
struct Params_map{ImgArray,ResultArray,N}
    # Specifies the dimensions 
    # by which to count including the negative direction
    rmask::Tuple{Vararg{Bool,N}}

    # size of result map
    result_size::Tuple{Vararg{Int,N}}

    # mirror buffers
    mirror_img::ImgArray
    mirror_result::ResultArray
end


"""
    Params_map(img, rmask)

Return parameters that can be used sequentially for all similar imgs.
"""
function Params_map(img::AbstractArray{T,N}, rmask::Tuple{Vararg{Bool,N}}) where T where N
    img_size = size(img)

    result_size = map((s, r) -> r ? 2s : s, img_size, rmask)

    mirror_img = similar(img)
    mirror_result = similar(img, Float64)

    Params_map(
        rmask,
        result_size,
        mirror_img,
        mirror_result
    )
end


function corr_function_map!(
    result::AbstractArray{T1,N},
    img::AbstractArray{T2,N}, 
    map_params::Params_map, 
    cf_params
) where T1 where T2 where N
    rmask = map_params.rmask
    mirror_img = map_params.mirror_img
    mirror_result = map_params.mirror_result

    for mask in BoolMask(N)
        if any(mask .& .!rmask)
            # skip because of symmetry
            continue
        end
        
        v_img = mirror(img, mask)
        v_result = maskres(result, img, mask, rmask)
        
        mirror_img .= v_img
        mirror_result .= 0
        
        correllation_function!(mirror_result, mirror_img, cf_params)
        v_result .= mirror_result
    end
    result
end


function corr_function_map(
    img::AbstractArray{T,N}, 
    CF_constructor; 
    parameters...
) where T where N
    cf_params, asymmetric = CF_constructor(img; parameters...)
    
    rmask = ones(Bool, N)
    rmask[N] = asymmetric
    rmask_tuple = Tuple(rmask)::Tuple{Vararg{Bool,N}}

    map_params = Params_map(img, rmask_tuple)
    result = similar(img, Float64, map_params.result_size)
    corr_function_map!(result, img, map_params, cf_params)
end
