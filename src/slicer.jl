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
