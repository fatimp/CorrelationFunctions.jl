using CorrelationFunctions
using XUnit

@testset "Diagonal slicing" begin include("diagonals.jl") end
@testset "Random data"      begin include("random.jl") end
@testset "Checkboard data"  begin include("checkboard.jl") end
@testset "Perlin noise"     begin include("perlin.jl") end
