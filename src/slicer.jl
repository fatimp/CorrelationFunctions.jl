# Where to put this function?
"""
    update_runs!(array, runs, n)

Add a sequence which start with the first element `runs` and decreases
by one to the first `n` elements of the vector `array`.
"""
function update_runs!(array :: Vector{Int}, runs, n)
    array[1:n] .+= Iterators.take(Iterators.countfrom(runs, -1), n)
end


function slice_generators(array :: Array{T,3},
                          direction :: Symbol) where T
    x, y, z = size(array)
    if direction == :x
        return (array[:,j,k] for j in 1:y for k in 1:z)
    elseif direction == :y
        return (array[i,:,k] for i in 1:x for k in 1:z)
    elseif direction == :z
        return (array[i,j,:] for i in 1:x for j in 1:y)
    else
        error("Unknown direction")
    end
end

with_doubling(iter, len) = imap(slice -> vcat(slice, slice[1:len]), iter)
