function crosscorr(a1, a2, direction, len, mode)
    @assert size(a1) == size(a2)
    check_direction(direction, a1, mode)
    result = zeros(Float64, len)
    plan = make_dft_plan(a1, mode, direction)

    slices1 = slices(a1, mode, direction)
    slices2 = slices(a2, mode, direction)
    for (slice1, slice2) in zip(slices1, slices2)
        padded1 = maybe_add_padding(slice1, mode)
        padded2 = maybe_add_padding(slice2, mode)

        ft1 = rfft_with_plan(padded1, plan)
        ft2 = rfft_with_plan(padded2, plan)
        ccft = @. ft1 * conj(ft2)
        cc = irfft_with_plan(ccft, length(padded1), plan)

        shifts = min(length(slice1), len)
        result[1:shifts] += cc[1:shifts]
    end

    return result
end

"""
    cross_correlation(array, phase1, phase2, direction[; len] [,mode = NonPeriodic()])

Calculate cross-correlation between `phase1` and `phase2` in
`array`. The meaning of optional arguments is the same as for `s2`
function.

See also: [`s2`](@ref), [`Utilities.AbstractDirection`](@ref),
[`Utilities.AbstractMode`](@ref).
"""
function cross_correlation(array, phase1, phase2, direction;
                           len  = (array |> size |> minimum) ÷ 2,
                           mode = NonPeriodic())
    ph1 = array .== phase1
    ph2 = array .== phase2

    m1 = maybe_apply_mask(ph1, mode)
    m2 = maybe_apply_mask(ph2, mode)
    cc = crosscorr(m1, m2, direction, len, mode)
    # m1 or m2 here, this does not matter
    return cc ./ normalization(m1, direction, len, mode)
end
