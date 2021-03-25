using CorrelationFunctions
using XUnit

@testset "Random data"     begin include("random.jl") end
@testset "Checkboard data" begin include("checkboard.jl") end
@testset "Perlin noise"    begin include("perlin.jl") end
