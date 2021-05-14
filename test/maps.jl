cut_direction(x :: Any, d :: Symbol) = cut_direction(x, Val(d))

cut_direction(array :: AbstractMatrix, :: Val{:x}) =
    let n = size(array, 2); array[n + 1:end, 1] end

cut_direction(array :: AbstractMatrix, :: Val{:y}) =
    let n = size(array, 2); array[n + 1, :] end

cut_direction(array :: AbstractMatrix, :: Val{:xy_main}) =
    let n = size(array, 2); [array[n + i, i] for i in 1:n] end

cut_direction(array :: AbstractMatrix, :: Val{:xy_anti}) =
    let n = size(array, 2); [array[n - i + 1, i] for i in 1:n] end

@testset "s2 map" begin
    testset = two_phase_noise_2d()
    for p in (false, true)
        directional = Directional.s2(testset, 1; len = 50, periodic = p,
                                     directions = [:x, :y, :xy_main, :xy_anti])
        cmap = Map.corr_function_map(testset, Map.s2; periodic = p)
        @test cut_direction(cmap, :x) ≈ directional[:x]
        @test cut_direction(cmap, :y) ≈ directional[:y]

        # TODO: Fix periodic diagonals
        if !p
            @test cut_direction(cmap, :xy_main) ≈ directional[:xy_main]
            @test cut_direction(cmap, :xy_anti) ≈ directional[:xy_anti]
        end
    end
end
