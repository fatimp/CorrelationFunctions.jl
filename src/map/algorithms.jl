function cnt_total(c :: AbstractArray)
    return map(size(c)) do s
        collect(s:-1:1)
    end
end

cnt_total_reshaped(c :: AbstractArray{T, 1}) where T = cnt_total(c)
cnt_total_reshaped(c :: AbstractArray{T, 2}) where T =
    let (t1, t2) = cnt_total(c);
        (reshape(t1, :, 1), reshape(t2, 1, :))
    end
cnt_total_reshaped(c :: AbstractArray{T, 3}) where T =
    let (t1, t2, t3) = cnt_total(c);
        (reshape(t1, :, 1, 1), reshape(t2, 1, :, 1), reshape(t3, 1, 1, :))
    end

function normalize_result(result   :: AbstractArray,
                          periodic :: Bool)
    local norm

    if periodic
        norm = result / length(result)
    else
        total = cnt_total_reshaped(result)
        norm = reduce(./, total; init = result)
    end

    return norm
end

function zeropad(image)
    s = size(image)
    s = (2 .* s) .- 1

    padded = similar(image, s)
    padded .= 0
    padded[axes(image)...] .= image
    return padded
end


function size_slice(img::AbstractArray{T, 3}, dim::Int) where T
    if dim == 1
        size(img)[2:3]
    elseif dim == 2
        (size(img, 1), size(img, 3))
    elseif dim == 3
        size(img)[1:2]
    else
        error("dim is not in {1, 2, 3}")
    end
end


function size_slice(img::AbstractArray{T, 2}, dim::Int) where T
    if dim == 1
        (size(img, 2),)
    elseif dim == 2
        (size(img, 1),)
    else
        error("dim is not in {1, 2}")
    end
end


function size_slice(img::AbstractVector, dim::Int)
    if dim == 1
        ()
    else
        error("dim is not 1")
    end
end
