function update_runs!(array, runs, n, :: NonPeriodic)
    array[1:n] .+= take(countfrom(runs, -1), n)
end

function update_runs!(array, runs, n, :: Periodic)
    array[1:n] .+= runs
end

# FIXME: For compatibility. Remove when masked computations are implemented
update_runs!(array, runs, n) = update_runs!(array, runs, n, NonPeriodic())

function update_runs_periodic!(array :: AbstractVector{T},
                               left  :: Integer,
                               right :: Integer,
                               len   :: Integer) where T <: Number
    sum = left + right
    nupdate = min(len, sum)
    f(n) = min(n - 1, left, right, sum - (n - 1))
    array[1:nupdate] .+= (f(n) for n in 1:nupdate)
end

function normalize_result(result, slices, mode :: AbstractMode)
    len = length(result)
    norm = zeros(Int, len)
    for slice in slices
        # Number of correlation lengths
        slen = length(slice)
        shifts = min(len, slen)

        update_runs!(norm, slen, shifts, mode)
    end

    return result ./ norm
end
