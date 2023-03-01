# Make a rotated cube
function make_cube(side, len)
    array = centered(falses(side, side, side))
    for idx in CartesianIndices(array)
        coord = 2 .* Tuple(idx) ./ side
        dist = coord .|> abs |> sum
        if dist < len
            array[idx] = true
        end
    end

    return array.parent
end

@testset "Check F_{sss} function for a cube" begin
    cube = make_cube(300, 0.45)
    @test_broken U.lowfreq_energy_ratio(square) > 0.97

    # F_{sss} for a cube is either undefined, zero or 2
    sss = D.surf3(cube, true; periodic = true, len = 20)
    slice = sss[D.PlaneXY()][15:end, 15:end] * 300^3
    @test maximum(relerr.(slice, 2)) < 0.2
end
