const known_directions = [U.DirX(),   U.DirY(),   U.DirZ(),
                          U.DirYZ(),  U.DirXZ(),  U.DirXY(),
                          U.DirZY(),  U.DirZX(),  U.DirYX(),
                          U.DirXYZ(), U.DirYXZ(), U.DirZYX()]
const axial_directions = [U.DirX(),   U.DirY(),   U.DirZ()]
const axial_directions_2d = [U.DirX(),   U.DirY()]

function mean_corrfn(fn, array, phase; directions = known_directions, kwargs...)
    cfs = map(dir -> fn(array, phase, dir; kwargs...), directions)
    return sum(cfs) / length(cfs)
end

"""
    draw_ball(s, r)

Draw an n-dimensional ball of radius R in an array of dimensionality
s.
"""
function draw_ball(s, r)
    array = zeros(Bool, s)
    center = @. (s ÷ 2) |> floor |> Int

    for idx in CartesianIndices(array)
        dist = sum((Tuple(idx) .- center) .^ 2)
        if dist <= r^2
            array[idx]	= true
        end
    end

    return array
end

isapprox_or_nan(x, y) = isapprox(x, y; atol = 1e-12) || (isnan(x) && isnan(y))

"""
    relerr(x, truex)

Calculate relative error between calculated value `x` and theoretical
value `truex`.
"""
relerr(x, truex) = abs((x - truex) / truex)

"""
    relerr_norm(x, truex)

Calculate |truex - x| / |truex|
"""
relerr_norm(x, truex) = norm(truex - x) / norm(truex)

"""
    scan(x)

Scan an array `x`, that is calculate running sum of its elements.
"""
function scan(x :: AbstractVector)
    result = copy(x)

    for i in 2:length(result)
        result[i] += result[i-1]
    end

    return result
end

"""
    integrate(f, range)

Integrate a function of one argument `f` across the range `range`
using Simpson's 1/3 rule
"""
function integrate(f :: Function, r :: AbstractRange)
    iter = zip(r, drop(r, 1))
    sum = 0.0

    for (a, b) in iter
        middle = (a + b) / 2
        sum += step(r)/6 * (f(a) + 4f(middle) + f(b))
    end

    return sum
end

# Value noise for tests

function lolrng(x    :: Integer,
                y    :: Integer,
                z    :: Integer,
                seed :: Integer)
    r1 = UInt32(x) * UInt32(0x1B873593)
    r2 = UInt32(y) * UInt32(0x19088711)
    r3 = UInt32(z) * UInt32(0xB2D05E13)

    r  = UInt32(seed) + r1 + r2 + r3
    r ⊻= r >> UInt32(5)
    r *= UInt32(0xCC9E2D51)
    return r / 0xffffffff
end

function octave(x    :: Number,
                y    :: Number,
                z    :: Number,
                oct  :: Integer,
                seed :: Integer)
    divisor = 2.0^(-oct)

    xi = (x / divisor) |> floor |> Int
    yi = (y / divisor) |> floor |> Int
    zi = (z / divisor) |> floor |> Int

    δx = rem(x, divisor) / divisor
    δy = rem(y, divisor) / divisor
    δz = rem(z, divisor) / divisor

    v000 = lolrng(xi,     yi,     zi, seed)
    v001 = lolrng(xi + 1, yi,     zi, seed)
    v010 = lolrng(xi,     yi + 1, zi, seed)
    v011 = lolrng(xi + 1, yi + 1, zi, seed)

    v100 = lolrng(xi,     yi,     zi + 1, seed)
    v101 = lolrng(xi + 1, yi,     zi + 1, seed)
    v110 = lolrng(xi,     yi + 1, zi + 1, seed)
    v111 = lolrng(xi + 1, yi + 1, zi + 1, seed)

    inter(v1, v2, x) = v1 + (v2 - v1)*x
    v00 = inter(v000, v001, δx)
    v01 = inter(v010, v011, δx)
    v10 = inter(v100, v101, δx)
    v11 = inter(v110, v111, δx)

    v0 = inter(v00, v01, δy)
    v1 = inter(v10, v11, δy)

    v = inter(v0, v1, δz)
    return v
end

value_noise(x        :: Number,
            y        :: Number,
            z        :: Number,
            octaves  :: Integer,
            seed     :: Integer) =
                mapreduce((x, octave) -> x ./ 2.0^octave, +,
                          (octave(x, y, z, o, seed) for o in 0:octaves-1),
                          countfrom(0)) ./ 2*(1 - 2.0^(-octaves))

two_phase_noise_3d() =
    let noise = [value_noise(x/10, y/10, z/10, 6, rand(UInt32))
                 for x in 1:50, y in 1:50, z in 1:50]
        noise .< 0.5
    end

two_phase_noise_2d() =
    let noise = [value_noise(x/10, y/10, 0, 6, rand(UInt32))
                 for x in 1:50, y in 1:50]
        noise .< 0.5
    end
