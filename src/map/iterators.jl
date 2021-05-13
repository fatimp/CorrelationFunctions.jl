# * BoolMask
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


# * Bresenham
struct Bresenham
    dx::Int
    dy::Int
    first_octant::Bool
end

Bresenham(x, y) = x >= y ? Bresenham(x, y, true) : Bresenham(y, x, false)

Base.length(B::Bresenham) = B.dx + 1
Base.eltype(::Type{Bresenham}) = Tuple{Int64,Int64}

function Base.iterate(B::Bresenham, state=(0, 0, 0))
    x, y, delta = state
    xnext, ynext = x, y
    
    if x > B.dx
        return
    end
    
    
    if x == 0
        # initialization of delta        
        delta = 2 * B.dy - B.dx
    else
        if delta < 0
            delta += 2 * B.dy
        else
            delta += 2 * B.dy - 2 * B.dx
        end
    end
    
    xnext += 1
    if delta >= 0
        ynext += 1
    end    
#     @show state
    
    if B.first_octant
        return ((x, y), (xnext, ynext, delta))
    else
        return ((y, x), (xnext, ynext, delta))
    end
end


# * NewBresenham
struct NewBresenham{N}
    dway::Tuple{Vararg{Int,N}}
    len::Int
end

NewBresenham(x::Tuple{Vararg{Int,N}}) where N = 
    NewBresenham{N}(x, maximum(abs.(x)))

Base.length(B::NewBresenham) = B.len + 1
Base.eltype(::Type{NewBresenham{N}}) where N = Tuple{Vararg{Int,N}}

function Base.iterate(B::NewBresenham, state::Int=0)
    if state > B.len
        return
    end
    
    past = state / B.len

    point = B.dway .* past
    
    (round.(Int, point), state + 1)
end
