function array_of_points(array :: AbstractArray{Bool})
    points = findall(array)
    return Array([idx[k] for k in 1:length(points[1]), idx in points]')
end

col_means(array :: AbstractMatrix) =
    Array(mean(array; dims = 1)')

function corrmatrix(array :: AbstractArray{Bool})
    points = array_of_points(array)
    m = col_means(points)
    dim = ndims(array)
    mat = mapreduce(+, axes(points, 1)) do idx
        coord = points[idx, :] - m
        map(Iterators.product(1:dim, 1:dim)) do (i, j)
            coord[i]*coord[j]
        end
    end

    return Symmetric(mat)
end


"""
    detect_anisotropy(array, phase)

Return a square matrix which characterizes anisotropy of the specified
phase in the input. Each column of the matrix is a unit vector. The
first vector points toward anisotropy and the other two vectors are
perpendicular to the first vector and to each other. The resulting
matrix can be used as an argument to `make_rotation`.

See also: [`rotate_array`](@ref), [`make_rotation`](@ref).
"""
function detect_anisotropy(array :: AbstractArray{<:Any, N}, phase) where N
    m = (array .== phase) |> corrmatrix |> eigvecs
    if det(m) < 0
        m[:,1] = -m[:,1]
    end

    return SMatrix{N, N}(m)
end
