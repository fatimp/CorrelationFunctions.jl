function cnt_total(img::AbstractArray, periodic::Bool)
    cnt_total(size(img), periodic)
end


function cnt_total(img::CuArray, periodic::Bool)
    cu(cnt_total(size(img), periodic))
end


function cnt_total(img_size::Tuple{Vararg{Int}}, periodic::Bool)
    if periodic
        result = *(img_size...)
    else
        result = zeros(Int, img_size)
        for ix in CartesianIndices(img_size)
            dir = Tuple(ix) .- 1
            result[ix] = *((img_size .- dir)...)
        end
    end
    result
end


function gradient_norm(img::Array)
    N = ndims(img)
    params = Tuple([true for _ in 1:N])
    sobels = [Kernel.sobel(params, i) for i in 1:N]

    dimgs = [imfilter(img, sobel) for sobel in sobels]

    map((x...) -> norm(x), dimgs...)
end


function gradient_norm(img::CuArray)
    imgCPU = Array(img)
    cu(gradient_norm(imgCPU))
end


"""
    cross_correlation!(f, g)

Compute cross correlation using FFT for xᵢ ≥ 0
and periodic boundary conditions.
Store result in `f`.
"""
function cross_correlation!(f, g)
    fft!(f)
    fft!(g)
    @. f = conj(f) * g
    ifft!(f)
end


"""
    self_correlation!(f)

Compute self correlation using FFT for xᵢ ≥ 0
and periodic boundary conditions
"""
function self_correlation!(f)
    fft!(f)
    @. f = abs2(f)
    ifft!(f)
end


function size_slice(img::AbstractArray{T, 3}, dim::Int) where T
    if dim == 1
        size(img)[2:3]
    elseif dim == 2
        (size(img, 1), size(img, 3))
    elseif dim == 3
        size(img)[1:2]
    else
        error("dim not in {1, 2, 3}")
    end
end


function size_slice(img::AbstractArray{T, 2}, dim::Int) where T
    if dim == 1
        (size(img, 2),)
    elseif dim == 2
        (size(img, 1),)
    else
        error("dim not in {1, 2}")
    end
end
