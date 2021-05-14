cut_direction(x :: Any, d :: Symbol) = cut_direction(x, Val(d))

cut_direction(array :: AbstractMatrix, :: Val{:x}) =
    let (x, y) = size(array); array[x÷2 + 1:end, 1] end

cut_direction(array :: AbstractMatrix, :: Val{:y}) =
    let (x, y) = size(array); array[x÷2 + 1, :] end

cut_direction(array :: AbstractMatrix, :: Val{:xy_main}) =
    let n = size(array, 2); [array[n + i, i] for i in 1:n] end

cut_direction(array :: AbstractMatrix, :: Val{:xy_anti}) =
    let n = size(array, 2); [array[n - i + 1, i] for i in 1:n] end

cut_direction(array :: AbstractArray{<:Any, 3}, :: Val{:x}) =
    let (x, y, z) = size(array); array[x÷2 + 1:end, y÷2 + 1, 1] end

cut_direction(array :: AbstractArray{<:Any, 3}, :: Val{:y}) =
    let (x, y, z) = size(array); array[x÷2 + 1, y÷2 + 1:end, 1] end

cut_direction(array :: AbstractArray{<:Any, 3}, :: Val{:z}) =
    let (x, y, z) = size(array); array[x÷2, y÷2, :] end

function testmap_2d(func, mapfunc)
    @testset "$func 2D map" begin
        testset = two_phase_noise_2d()
        for p in (false, true)
            # TODO: Fix periodic diagonals
            directions = p ? [:x, :y] : [:x, :y, :xy_main, :xy_anti]
            directional = func(testset, 1; len = 50, periodic = p, directions = directions)
            cmap = mapfunc(testset, p)

            for (direction, data) in directional
                # ! FIXME l2 on diagonals are broken
                if func == Directional.l2 && direction ∈ [:xy_main, :xy_anti] 
                    @test_broken cut_direction(cmap, direction) ≈ data
                else
                    @test cut_direction(cmap, direction) ≈ data
                end
            end
        end
    end
end

mapfn(map) = (x, p) -> Map.corr_function_map(x, map; periodic = p)
testmap_2d(Directional.l2,       mapfn(Map.l2))
testmap_2d(Directional.s2,       mapfn(Map.s2))
testmap_2d(Directional.c2,       mapfn(Map.c2))
testmap_2d(Directional.surfsurf, mapfn(Map.ss))
testmap_2d(Directional.surfvoid,
           (x -> let n = size(x, 2); x[:, n÷2 + 1:end] end) ∘ mapfn(Map.sv))

function testmap_3d(func, mapfunc)
    @testset "$func 3D map" begin
        testset = two_phase_noise_3d()
        for p in (false, true)
            directional = func(testset, 1; len = 50, periodic = p)
            cmap = mapfunc(testset, p)

            for (direction, data) in directional
                @test cut_direction(cmap, direction) ≈ data
            end
        end
    end
end

# Takes a whole eternity to execute
#testmap_3d(Directional.l2,       mapfn(Map.l2))
testmap_3d(Directional.s2,       mapfn(Map.s2))
# Takes a whole eternity to execute
#testmap_3d(Directional.c2,       mapfn(Map.c2))
testmap_3d(Directional.surfsurf, mapfn(Map.ss))
testmap_3d(Directional.surfvoid,
           (x -> let n = size(x, 3); x[:, :, n÷2 + 1:end] end) ∘ mapfn(Map.sv))
