using CorrelationFunctions
using XUnit
using PoissonRandom

@testset "Diagonal slicing"        begin include("diagonals.jl") end
@testset "Random data"             begin include("random.jl") end
@testset "Checkboard data"         begin include("checkboard.jl") end
@testset "Perlin noise"            begin include("perlin.jl") end
@testset "Non-overlapping spheres" begin include("poisson-spheres.jl") end
