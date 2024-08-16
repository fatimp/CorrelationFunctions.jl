map_slice_len(dim :: Int64, :: Periodic) = dim
map_slice_len(dim :: Int64, :: NonPeriodic) = dim ÷ 2 + 1

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
    # This works only for arrays with all dimensions being odd, which
    # is OK becase we usually pad arrays to dimensions 2s[i]-1 for
    # non-periodic computations
    @assert all(isodd, size(c))

    # cnt_total[i,j,k] = cnt_total_axis_1[i] * cnt_total_axis_2[j] * …
    result = mapreduce(.*, axes(c), size(c), 1:N) do ix, es, dim
        shape = collect(i == dim ? (:) : 1 for i in 1:N)
        s = (es + 1) ÷ 2
        # We get a descending sequence s, s-1, s-2, …, 1 followed by
        # an ascending sequence 1, 2, …, s-1
        cnt = @. ifelse(ix <= s, s - ix + 1, ix - s)
        # It's necessary to give a hint what a result of reshape is
        reshape(cnt, shape...) :: Array{Int64, N}
    end

    return maybe_upload_to_gpu(result, c)
end

normalize_result(result, :: Periodic)    = result  / length(result)
normalize_result(result, :: NonPeriodic) = result ./ cnt_total(result)

function normalize_result(result, mode :: Mask)
    mask = mode.mask
    n = autocorr(mask, NonPeriodic())

    return result ./ n
end

maybe_apply_mask(array, :: AbstractMode) = array
function maybe_apply_mask(array, mode :: Mask)
    mask = mode.mask
    @assert size(array) == size(mask)
    return array .* mask
end
