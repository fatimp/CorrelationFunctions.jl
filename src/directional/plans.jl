abstract type AbstractDFTPlan end

struct BogusDFTPlan <: AbstractDFTPlan end

struct FFTWPlan{F, I} <: AbstractDFTPlan
    forward :: F
    inverse :: I
end

# XXX: Remove this function after refactoring of slicer.jl
topology_to_periodic(:: Torus) = true
topology_to_periodic(:: Plane) = false

function make_dft_plan(array     :: AbstractArray,
                       topology  :: AbstractTopology,
                       direction :: AbstractDirection)
    periodic = topology_to_periodic(topology)
    check_direction(direction, array, periodic)

    same_length = slices_have_same_length(topology, direction)
    if same_length
        slice = similar(first(slice_generators(array, periodic, direction)), Int)
        padded = maybe_pad_with_zeros(slice, topology)
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
