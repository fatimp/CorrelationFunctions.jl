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

function normalization(array, direction, len, mode :: AbstractMode)
    norm = zeros(Int, len)
    for slice in slices(array, mode, direction)
        # Number of correlation lengths
        slen = length(slice)
        shifts = min(len, slen)

        update_runs!(norm, slen, shifts, mode)
    end

    return norm
end

normalization(array, direction, len, mode :: Mask) =
    autocorr(mode.mask, direction, len, NonPeriodic())
