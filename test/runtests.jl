using CorrelationFunctions
using Test

# Test correlation functions averaged on all directions
ms2 = mean ∘ s2
ml2 = mean ∘ l2

# The same with periodic boundary conditions
pms2 = mean ∘ (array, len, phase) -> s2(array, len, phase; periodic = true)
pml2 = mean ∘ (array, len, phase) -> l2(array, len, phase; periodic = true)

@testset "Random data"     begin include("random.jl") end
@testset "Checkboard data" begin include("checkboard.jl") end
