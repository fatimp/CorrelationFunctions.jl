"""
    maybe_upload_to_gpu(a1, a2)

Construct `CuArray` from `a1` if `a2` is `CuArray`.
"""
function maybe_upload_to_gpu end

maybe_upload_to_gpu(array :: AbstractArray, :: AbstractArray) = array
maybe_upload_to_gpu(array :: AbstractArray, :: CuArray) = CuArray(array)
