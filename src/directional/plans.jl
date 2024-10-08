abstract type AbstractDFTPlan end

struct BogusDFTPlan <: AbstractDFTPlan end

struct FFTWPlan{F, I} <: AbstractDFTPlan
    forward :: F
    inverse :: I
end

function make_dft_plan(array     :: AbstractArray,
                       mode      :: AbstractMode,
                       direction :: AbstractDirection)
    check_direction(direction, array, mode)

    same_length = slices_have_same_length(mode, direction)
    if same_length
        slice = similar(first(slices(array, mode, direction)), Int)
        padded = maybe_add_padding(slice, mode)
        len = length(padded)
        fwd = plan_rfft(padded)
        ft  = fwd * padded
        inv = plan_irfft(ft, len)
        return FFTWPlan(fwd, inv)
    else
        return BogusDFTPlan()
    end
end

rfft_with_plan(array :: AbstractVector, :: BogusDFTPlan) = rfft(array)
rfft_with_plan(array :: AbstractVector, plan :: FFTWPlan) = plan.forward * array

irfft_with_plan(array :: AbstractVector, len :: Integer, :: BogusDFTPlan) = irfft(array, len)
irfft_with_plan(array :: AbstractVector, ::Any, plan :: FFTWPlan) = plan.inverse * array
