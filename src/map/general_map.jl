function mirror(image, mask)
    ixs = map((ix, m) -> m ? reverse(ix) : ix, axes(image), mask)
    view(image, ixs...)
end

function corr_function_map(img :: AbstractArray, cf_params)
    result = CFMap(img)

    mirror_img = similar(img)
    mirror_result = similar(img, Float32)

    for mask in BoolMask(ndims(img))
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
