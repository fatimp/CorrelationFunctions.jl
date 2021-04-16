# CorrelationFunctions.jl
[![CI](https://github.com/shamazmazum/CorrelationFunctions.jl/actions/workflows/test.yml/badge.svg)](https://github.com/shamazmazum/CorrelationFunctions.jl/actions/workflows/test.yml)

Short how-to in examples:

~~~~{.jl}
using CorrelationFunctions

# Generate random 200x250x150 array of zeros and ones
a = rand(0:1, (200, 250, 150));
~~~~

Now you can calculate, for example, `L2(x)` for the phase `1` of the array `a`
where `x` runs from `1` to `10`. `len` is an optional argument and defaults to
half of the minimal direction of the array (`75` in this case).

~~~~{.jl}
l = l2(a, 1; len = 10)
┌────────────┬─────────────┬─────────────┐
│          x │           y │           z │
├────────────┼─────────────┼─────────────┤
│   0.500067 │    0.500067 │    0.500067 │
│   0.250169 │    0.250069 │    0.249895 │
│   0.125177 │    0.125132 │    0.124869 │
│   0.062549 │   0.0626169 │   0.0624302 │
│  0.0312185 │   0.0313031 │   0.0312011 │
│  0.0155847 │   0.0156507 │   0.0155863 │
│ 0.00777526 │  0.00784413 │  0.00777139 │
│ 0.00388573 │  0.00393841 │  0.00388657 │
│ 0.00193833 │  0.00198802 │  0.00194197 │
│  0.0009526 │ 0.000998479 │ 0.000962695 │
└────────────┴─────────────┴─────────────┘
~~~~

You can average the calculated data along multiple
directions:

~~~~{.jl}
mean(l, [:x, :y])
10-element Array{Float64,1}:
 0.5000665333333333
 0.250118935208438
 0.12515452404978136
 0.06258300388579152
 0.03126089613034623
 0.015617800511508951
 0.007809797875984926
 0.003912166580622957
 0.001963278008298755
 0.0009756470383880493
~~~~

Calling `mean` without the second argument averages data along all directions.

By default correlation functions are calculated along array axes. You can
calculate them along diagonals too, for example, let's calculate `L2` along
`(1,1,0)` and `(1, -1, 0)` directions.

~~~~{.jl}
l = l2(a, 1; len = 10, directions = [:xy_main, :xy_anti])
┌─────────────┬─────────────┐
│     xy_main │     xy_anti │
├─────────────┼─────────────┤
│    0.500067 │    0.500067 │
│    0.250199 │    0.250186 │
│    0.125144 │    0.125092 │
│   0.0626356 │   0.0625594 │
│   0.0313083 │   0.0313232 │
│   0.0156661 │   0.0156908 │
│  0.00779083 │  0.00785716 │
│  0.00385538 │   0.0039246 │
│   0.0019037 │  0.00195348 │
│ 0.000946319 │ 0.000967609 │
└─────────────┴─────────────┘
~~~~

The result is basically the same, because zeros and ones are distributed
uniformly in the array.

You can impose periodic boundary conditions onto `array` by calling `L2` with
`periodic = true` argument.

~~~~{.jl}
l = CorrelationFunctions.l2(a, 1; len = 10, periodic = true)
┌─────────────┬─────────────┬─────────────┐
│           x │           y │           z │
├─────────────┼─────────────┼─────────────┤
│    0.500043 │    0.500087 │     0.49997 │
│     0.25015 │    0.250094 │    0.249827 │
│     0.12515 │    0.125151 │    0.124822 │
│   0.0625262 │   0.0626342 │   0.0623815 │
│   0.0311995 │   0.0313199 │    0.031155 │
│   0.0155828 │   0.0156675 │   0.0155542 │
│   0.0077851 │  0.00786115 │  0.00774896 │
│  0.00389583 │  0.00394862 │  0.00387399 │
│  0.00194812 │  0.00198981 │  0.00192987 │
│ 0.000958806 │ 0.000999336 │ 0.000954834 │
└─────────────┴─────────────┴─────────────┘
~~~~

The other two functions, `s2` and `c2` have the same interface, so I will not
describe them in details.
