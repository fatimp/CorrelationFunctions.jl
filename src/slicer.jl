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

struct PeriodicIterator
    iterator
    len :: Integer
end

add_tail(:: Int, :: Nothing) = nothing
function add_tail(len :: Int, it :: Tuple{Vector, Any})
    slice, state = it
    return vcat(slice, slice[1:len-1]), state 
end

Base.iterate(iter :: PeriodicIterator) =
    add_tail(iter.len, iterate(iter.iterator))
Base.iterate(iter :: PeriodicIterator, state ) =
    add_tail(iter.len, iterate(iter.iterator, state))
