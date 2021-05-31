function mirror(image, mask)
    ixs = map((ix, m) -> m ? reverse(ix) : ix, axes(image), mask)
    view(image, ixs...)
end


"""
    Params_map

Store parameters and inner pre-allocated matrices.
"""
struct Params_map{ImgArray,ResultArray}
    mirror_img::ImgArray
    mirror_result::ResultArray
end


"""
    Params_map(img, rmask)

Return parameters that can be used sequentially for all similar imgs.
"""
function Params_map(img)
    mirror_img = similar(img)
    mirror_result = similar(img, Float64)

    Params_map(
        mirror_img,
        mirror_result
    )
end


function corr_function_map!(
    result::CFMap{T1, N},
    img::AbstractArray{T2,N},
    map_params::Params_map,
    cf_params
) where T1 where T2 where N
    mirror_img = map_params.mirror_img
    mirror_result = map_params.mirror_result

    for mask in BoolMask(N)
        if result.cf_type == :central_symmetry && mask[end]
            continue
        elseif result.cf_type == :periodic_point_point && any(mask)
            continue
        end


        mirror_img .= mirror(img, mask)
        mirror_result .= 0

        correllation_function!(mirror_result, mirror_img, cf_params)

        v_result = cut_cfmap(result, mask)
        v_result .= mirror_result
    end
    result
end


"""
    corr_function_map(img, CF_constructor; parameters...)

Compute correlation function map in all meaningfull directions.

CF_constructor ∈ {Map.s2, Map.l2, Map.c2, Map.ss, Map.sv}
"""
function corr_function_map(
    img::AbstractArray{T,N},
    CF_constructor;
    parameters...
) where T where N
    cf_params, cf_type = CF_constructor(img; parameters...)

    map_params = Params_map(img)
    result = CFMap(img, cf_type)

    corr_function_map!(result, img, map_params, cf_params)
end
