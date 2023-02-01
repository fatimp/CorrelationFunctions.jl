function s3 end

function s3_at_point(array :: AbstractArray,
                     x     :: AbstractVector,
                     y     :: AbstractVector)
    shift1 = circshift(array, x)
    shift2 = circshift(array, y)
    mul = @. array * shift1 * shift2
    return sum(mul)
end

function s3_plane(array :: AbstractArray,
                  plane :: AbstractPlane,
                  len)
    shift1, shift2 = unit_shifts(array, plane)
    result = zeros(Int, (len, len))

    for idx in CartesianIndices(result)
        s1 = (idx[1] - 1) * shift1
        s2 = (idx[2] - 1) * shift2
        result[idx] = s3_at_point(array, s1, s2)
    end

    return result / length(array)
end

function s3(array  :: AbstractArray;
            planes :: Vector{AbstractPlane} = default_planes(array),
            len                             = (array |> size |> minimum) ÷ 2)
    calc_s3 = plane -> plane => s3_plane(array, plane, len)
    return Dict{AbstractPlane, Matrix{Float64}}(map(calc_s3, planes))
end
