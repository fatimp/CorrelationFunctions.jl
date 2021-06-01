"""
    BoolMask(n::Int)

Iterate over Boolean representaion of 1:2^n - 1
"""
struct BoolMask
    N::Int
end

Base.length(BM::BoolMask) = 2^BM.N
Base.eltype(::Type{BoolMask}) = Array{Bool,1}

function Base.iterate(BM::BoolMask, state::Int=0)
    n = BM.N
    if state == 2^n
        return
    end

    mask = digits(Bool, state, base=2, pad=n)

    return (mask, state + 1)
end
