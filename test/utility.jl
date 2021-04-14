"""
    relerr(x, truex)

Calculate relative error between calculated value `x` and theoretical
value `truex`.
"""
relerr(x, truex) = abs((x - truex) / truex)

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
