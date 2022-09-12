function mirror(image, mask)
    ixs = map((ix, m) -> m ? reverse(ix) : ix, axes(image), mask)
    view(image, ixs...)
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
    cf_params = CF_constructor(img; parameters...)
    result = CFMap(img)

    mirror_img = similar(img)
    mirror_result = similar(img, Float32)

    for mask in BoolMask(N)
        if mask[end]
            continue
        end

        mirror_img .= mirror(img, mask)
        mirror_result .= 0

        correllation_function!(mirror_result, mirror_img, cf_params)

        v_result = cut_cfmap(result, mask)
        v_result .= mirror_result
    end

    return result
end
