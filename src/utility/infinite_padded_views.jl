# This is a wrapper around arrays which provides "infinite" PaddedView
# interface: in-bounds access to an array works as usual and
# out-of-bounds access returns a special "fill" value.

struct InfinitePaddedView{T,N,A} <: AbstractArray{T,N}
    parent :: A
    fill   :: T
end

InfinitePaddedView(array :: A, fill) where A <: AbstractArray =
    InfinitePaddedView{eltype(A), ndims(A), A}(array, eltype(A)(fill))

Base.size(array :: InfinitePaddedView) = size(array.parent)

function Base.getindex(array :: InfinitePaddedView{<: Any, N},
                       idx   :: Vararg{Int, N}) where N
    if checkbounds(Bool, array.parent, idx...)
        return array.parent[idx...]
    else
        return array.fill
    end
end

function Base.setindex!(array :: InfinitePaddedView{<: Any, N},
                        x,
                        idx   :: Vararg{Int, N}) where N
    array.parent[idx...] = x
end

Base.axes(array :: InfinitePaddedView) = axes(array.parent)

Base.similar(array :: InfinitePaddedView) =
    InfinitePaddedView(similar(array.parent), array.fill)
