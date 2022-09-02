"""
    dir_from_map(m, dir)

Extract a direction from a correlation map. `direction` can be either
`:x`, `:y`, `:z`, `:xy` or `:yx`.
"""
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

"""
    average_directions(cfmap)

Average correlation map `cfmap` over all directions. The result is a
vector with indices being equal to correlation length + 1.
"""
function average_directions(cfmap :: AbstractArray{T}) where T
    num  = Dict{Int, Int}()
    vals = Dict{Int, Vector{T}}()
    for idx in CartesianIndices(cfmap)
        dist = ((Tuple(idx) .- 1) .^ 2) |> sum |> sqrt |> round |> Int
        v = get(vals, dist, T[])
        push!(v, cfmap[idx])
        num[dist]  = get(num, dist, 0) + 1
        vals[dist] = v
    end

    len = num |> keys |> maximum
    avg = Vector{Float64}(undef, len)

    for idx in 1:len
        n = get(num,  idx-1, 1)
        v = get(vals, idx-1, T[])
        avg[idx] = sum(v) / n
    end

    return avg
end

function cnt_total_(c; periodic=false, original=false)
    if periodic
        [length(c)]
    else
        if original
            ixes = map(size(c)) do s
                collect(s:-1:1)
            end
        else
            ixes = map(axes(c), size(c)) do ix, es
                s = (es + 1) ÷ 2
                cnt = @. s - abs(ix - s)
                ifftshift(cnt)
            end
        end
        map(ixes, 1:length(ixes)) do cnt, d
            if d == 1
                cnt
            elseif d == 2
                reshape(cnt, 1, :)
            elseif d == 3
                reshape(cnt, 1, 1, :)
            end
        end
    end
end


cnt_total(c; periodic=false, original=false) = cnt_total_(c; periodic, original)
cnt_total(c::CuArray; periodic=false, original=false) = cnt_total_(c; periodic, original) .|> cu


function zeropad(image)
    s = size(image)
    s = (2 .* s) .- 1

    padded = similar(image, s)
    padded .= 0
    padded[axes(image)...] .= image
    return padded
end
