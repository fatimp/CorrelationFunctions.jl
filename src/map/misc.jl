map_slice_len(dim :: Int64, :: Periodic) = dim
map_slice_len(dim :: Int64, :: AbstractMode) = dim ÷ 2

"""
    dir_from_map(m, dir)

Extract a direction from a correlation map. `direction` must be of
type `AbstractDirection`.

See also: [`Utilities.AbstractDirection`](@ref).
"""
function dir_from_map end

function dir_from_map(x    :: AbstractArray, :: DirX;
                      mode :: AbstractMode = NonPeriodic())
    q = map_slice_len(size(x, 1), mode)
    return x[CartesianIndex.(1:q, 1, 1)]
end

function dir_from_map(x    :: AbstractArray, :: DirY;
                      mode :: AbstractMode = NonPeriodic())
    q = map_slice_len(size(x, 2), mode)
    return x[CartesianIndex.(1, 1:q, 1)]
end

function dir_from_map(x    :: AbstractArray, :: DirZ;
                      mode :: AbstractMode = NonPeriodic())
    q = map_slice_len(size(x, 3), mode)
    return x[CartesianIndex.(1, 1, 1:q)]
end

function dir_from_map(x    :: AbstractArray, :: DirXY;
                      mode :: AbstractMode = NonPeriodic())
    q = map_slice_len(minimum(size(x)[1:2]), mode)
    return x[CartesianIndex.(1:q, 1:q, 1)]
end

function dir_from_map(x    :: AbstractArray, :: DirYX;
                      mode :: AbstractMode = NonPeriodic())
    a = minimum(size(x)[1:2])
    q = map_slice_len(a, mode)
    b = 1:q
    c = @. mod(1 - b, a) + 1
    return x[CartesianIndex.(c, b, 1)]
end

"""
    average_directions(cfmap; len = (cfmap |> size |> minimum) ÷ 2)

Average correlation map `cfmap` over all directions. The result is a
vector of length `len` with indices being equal to correlation length + 1.
"""
function average_directions(cfmap :: AbstractArray{T};
                            len   :: Integer = (cfmap |> size |> minimum) ÷ 2) where T
    counter = zeros(Int, len)
    accum   = zeros(T,   len)

    for idx in CartesianIndices(cfmap)
        dist = round(Int, ((Tuple(idx) .- 1) .^ 2) |> sum |> sqrt) + 1
        if dist <= len
            accum[dist] += cfmap[idx]
            counter[dist] += 1
        end
    end

    return accum ./ counter
end

function cnt_total(c :: AbstractArray{<:Any, N}) where N
    # This works only for arrays with all dimensions being even, which
    # is OK becase we pad arrays to dimensions 2s[i] for non-periodic
    # computations.
    @assert all(iseven, size(c))

    # cnt_total[i,j,k] = cnt_total_axis_1[i] * cnt_total_axis_2[j] * …
    result = mapreduce(.*, axes(c), size(c), 1:N) do ix, es, dim
        shape = collect(i == dim ? (:) : 1 for i in 1:N)
        s = es ÷ 2
        # We get a descending sequence s, s-1, s-2, …, 0 followed by
        # an ascending sequence 1, 2, …, s-1
        cnt = @. ifelse(ix <= s + 1, s - ix + 1, ix - s - 1)
        # It's necessary to give a hint what a result of reshape is
        reshape(cnt, shape...) :: Array{Int64, N}
    end

    return maybe_upload_to_gpu(result, c)
end

function autocorr(array)
    ft = rfft(array)
    # There is no method irfft(:: CuArray{Float64}, :: T)!
    # Do not use abs2 here or in Map.c2!
    return irfft(ft .* conj.(ft), size(array, 1))
end

function crosscorr(a1, a2)
    @assert size(a1) == size(a2)
    s = size(a1, 1)
    plan = plan_rfft(a1)

    ft1 = plan * a1
    ft2 = plan * a2
    ccf = @. ft1 * conj(ft2)
    return irfft(ccf, s)
end

normalize_result(result, :: Periodic)    = result  / length(result)
function normalize_result(result, :: NonPeriodic)
    na = cnt_total(result)
    return @. ifelse(na == 0, NaN, result / na)
end

function normalize_result(result, mode :: Mask)
    padded = maybe_add_padding(mode.mask, NonPeriodic())
    ac = round.(autocorr(padded))
    n = @. ifelse(ac == 0, NaN, ac)
    return result ./ n
end
