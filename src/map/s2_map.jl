struct Params_S2{ComplexArray,Total}
    # boundary conditions
    periodic::Bool
    # normalization
    total::Total

    # algorithm-specific

    # fft buffers
    complex_img::ComplexArray
end


function Params_S2(img; periodic::Bool=true)
    box = size(img)
    complex_box = periodic ? box : box .* 2

    total = cnt_total(img, periodic)

    p = Params_S2(
        periodic,
        total,
        similar(img, ComplexF32, complex_box),
    )
    cf_type = periodic ? :periodic_point_point : :central_symmetry
    p, cf_type
end


"""
    correllation_function!(res, img, params::Params_S2)

Compute S₂ correlation function in positive directions
"""
function correllation_function!(res, img, params::Params_S2)
    ix = CartesianIndices(img)

    f = params.complex_img .= 0

    v_f = view(f, ix)
    v_f .= img

    self_correlation!(f)

    res .= real.(v_f) ./ params.total
end


"""
    s2(image; periodic = false)

Calculate `S₂` (two point) correlation function map
for the N-dimensional image and return a `CFMap` object.

The `image` contains the probability of the voxel being in the correct phase.

# Examples
```jldoctest
julia> s2([1 0; 0 1]; periodic=true).result
2×2 Matrix{Float32}:
 0.5  0.0
 0.0  0.5
```
"""
function s2(image; periodic::Bool=false)
    corr_function_map(image, Params_S2; periodic)
end


function s2(image, phase::Int; periodic::Bool=false)
    one_phase_img = image .== phase
    corr_function_map(one_phase_img, Params_S2; periodic)
end


function cross_correlation(a, b = a)
    A = rfft(a)
    if b === a
        C = @. abs2(A)
    else
        B = rfft(b)
        C = @. A * conj(B)
    end
    irfft(C, size(a, 1))
end


function cnt_total(c; periodic=false)
    if periodic
        length(c)
    else
        orig_size = [(s + 1) ÷ 2 for s in size(c)]
        ixes = map((ix, s) -> ix .- s, axes(c), orig_size)
        qshift = map(Iterators.product(ixes...)) do ix
            mapreduce(*, orig_size, ix) do s, i
                s - abs(i)
            end
        end
        ifftshift(qshift)
    end
end


function expand(image)
    esize = [2s - 1 for s in size(image)]
    e_image = similar(image, esize...) .= false
    e_image[axes(image)...] .= image
    e_image
end


"""
    s2(image; periodic = false)

Calculate `S₂` (two point) correlation function map
for the N-dimensional image and return a `CFMap` object.

The `image` contains the probability of the voxel being in the correct phase.

# Examples
```jldoctest
julia> s2([1 0; 0 1]; periodic=true).result
2×2 Matrix{Float32}:
 0.5  0.0
 0.0  0.5
```
"""
function new_s2(image; periodic=false)
    if periodic
        A = image
    else
        A = expand(image)
    end
    c = cross_correlation(A)
    q = cnt_total(c; periodic)
    c ./ q
end
