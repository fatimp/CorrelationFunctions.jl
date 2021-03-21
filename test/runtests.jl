using CorrelationFunctions
using Test

# Test correlation functions averaged on all directions
ms2 = mean ∘ s2
ml2 = mean ∘ l2

@testset "Random data"     begin include("random.jl") end
@testset "Checkboard data" begin include("checkboard.jl") end
