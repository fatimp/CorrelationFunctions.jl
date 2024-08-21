"""
    update_runs!(array, runs, n)

Add a sequence which start with the first element `runs` and decreases
by `dec` to the first `n` elements of the vector `array`.
"""
function update_runs!(array :: AbstractVector{T}, runs, n) where T <: Number
    array[1:n] .+= take(countfrom(runs, -1), n)
end

function update_runs_periodic!(array :: AbstractVector{T},
                               left  :: Integer,
                               right :: Integer,
                               len   :: Integer) where T <: Number
    sum = left + right
    nupdate = min(len, sum)
    f(n) = min(n - 1, left, right, sum - (n - 1))
    array[1:nupdate] .+= (f(n) for n in 1:nupdate)
end
