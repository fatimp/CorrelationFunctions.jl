# Accessing Data

The most functions in this package (with exception to `pore_size` and
`chord_length`) return a value of type `CorrelationData`:

```@example
using CorrelationFunctions
using Random

a = l2(rand(MersenneTwister(1453), 0:1, (100, 100, 100)), 1)
```

You can extract the values along any computed direction using indexing operator:
```@example
using CorrelationFunctions
using Random

a = l2(rand(MersenneTwister(1453), 0:1, (100, 100, 100)), 1)
a[:y]
```

Also you can average results along multiple directions using `StatsBase.mean`
function:
```@example
using CorrelationFunctions
using Random
using StatsBase

a = l2(rand(MersenneTwister(1453), 0:1, (100, 100, 100)), 1)
mean(a, [:x, :y])
```

Calling `StatsBase.mean` without the second argument averages along all computed
directions.
