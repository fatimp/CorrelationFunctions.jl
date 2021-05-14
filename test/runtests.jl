using CorrelationFunctions
using XUnit
using PoissonRandom
using Base.Iterators: zip, drop, countfrom
using StatsBase: mean

# Test utilities
include("utility.jl")

@testset "Diagonal slicing"  begin include("diagonals.jl") end
@testset "Short sequences"   begin include("short.jl") end
@testset "Random data"       begin include("random.jl") end
@testset "Checkboard data"   begin include("checkboard.jl") end
@testset "Value noise"       begin include("value-noise.jl") end
@testset "Overlapping disks" begin include("poisson-spheres-2d.jl") end
@testset "Overlapping balls" begin include("poisson-spheres-3d.jl") end
