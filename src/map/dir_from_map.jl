function dir_from_map(m::AbstractArray, direction; periodic=false)
    if direction == :x
        q = periodic ? size(m, 1) : size(m, 1) ÷ 2 + 1
        ixs = CartesianIndex.(1:q, 1, 1)
    elseif direction == :y
        q = periodic ? size(m, 2) : size(m, 2) ÷ 2 + 1
        ixs = CartesianIndex.(1, 1:q, 1)
    elseif direction == :z
        q = periodic ? size(m, 3) : size(m, 3) ÷ 2 + 1
        ixs = CartesianIndex.(1, 1, 1:q)
    elseif direction == :xy
        a = minimum(size(m)[1:2])
        q = periodic ? a : a ÷ 2 + 1
        ixs = CartesianIndex.(1:q, 1:q, 1)
    elseif direction == :yx
        a = minimum(size(m)[1:2])
        q = periodic ? a : a ÷ 2 + 1
        b = 1:q
        c = @. mod(1 - b, a) + 1
        ixs = CartesianIndex.(c, b, 1)
    end
    m[ixs]
end
