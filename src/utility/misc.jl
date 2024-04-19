"""
    maybe_upload_to_gpu(a1, a2)

Construct `CuArray` from `a1` if `a2` is `CuArray`.
"""
function maybe_upload_to_gpu end

maybe_upload_to_gpu(array :: AbstractArray, :: AbstractArray) = array
maybe_upload_to_gpu(array :: AbstractArray, :: CuArray) = CuArray(array)

"""
    check_rank(a, dim)

Check if an array `a` has dimensionality equal or greater than
`dim`. Return `a` or signal an error.
"""
function check_rank(a, dim)
    if ndims(a) ≥ dim
        return a
    else
        error("This function is only defined for arrays of dimensionality ≥ $(dim)")
    end
end
