# Changelog

## Version 0.14.0

* Enhancement: Start to use `label_components` from ImageMorphology.jl again
  when it does not crash Julia.
* Incompatible change: When computing correlation functions in non-periodic
  mode, pad to 2*dim[k] instead of 2*dim[k]-1. This makes computational speed of
  Fourier transform more predictable. Meaningless elements of correlation maps
  will contain `NaN`s.

## Version 0.13.0

Incompatible changes:

* `make_pattern()` is renamed to `right_triangles()` (as it creates a
  pattern based on right triangles).
* `AbstractIndicator` type and all its subtypes are removed.
* Image rotation stuff is completely removed.
* `periodic` keyword argument is replaced with `mode`. A mode is a
  value of type `AbstractMode`. Computation using a mask can be
  performed by specifying `mode = Mask(mask)`.

Improvements:

* Allow arument to three-point functions to be negative.

Bug fixes:

* Fixed normalization of three-point functions with arbitrary pattern.

## Version 0.12.0

* Incompatible change: all directional 2-point correlation functions
  now require a direction as a mandatory third argument and return
  an array. CorrelationData structure is removed.
* Incompatible change: S2FTPlans structure along with `plans` argument
  are removed.
* Improvement: 3-point CFs can be calculated in arbitrary points.
  They accept an array of points as the third and the fourth arguments.
* Optimization: Directional.s2 is faster when calculated for diagonal
  directions in non-periodic mode.

## Version 0.11.1

* Optimization: Import time is reduced by 30%
* Incompatible change: Surface correlation functions are fool-proofed
  against input arrays of incompatible dimensionality.

## Version 0.11.0

Incompatible changes:

* `Directional.chord_length` and `Directional.pore_size` now return
  arrays of chord_length and pore sizes rather than a histogram.

## Version 0.10.3

New features:

* Add a new function `Utilities.rotate_array` which performs image
  rotation. Like `imrotate`, but works also for 3D images. It has
  support for coordinate transforms in a matrix form and rotations
  around a vector.
* Add a new function `Utilities.detect_anisotropy` which returns
  a distinctive direction for a N-dimensional image, as well as
  another N-1 vectors which form a basic. This basic may be used
  for the previous function.

Removed functions and methods:

* `Utilities.read_cuboid` method which required JSON metadata was
removed. This package does not depend on `JSON` now.

Other stuff:

* Documentation actualization.

## Version 0.10.2

????

## Version 0.10.1

* Remove erroneous dependency on Documenter.jl
* Faster Directional.c2 for periodic mode. Also this function consumes
  less memory.
* Code style fixes.

## Version 0.10.0

?????

## Version 0.9.0

Changes from the previous version:

* Symbolic designators for directions are replaced with special
  types (subtypes of AbstractDirection).
* CorrelationData implements AbstractDict interface now and is
  a dictionary with keys of type AbstractDirection and values of
  type Vector{Float64}.
* Three-point correlation functions `s3` and `c3` (cluster function)
  are added.

## Version 0.8.0

* Function `Map.average_directions` is much faster now.
* Function `Utilities.extract_edges` uses `imfilter` from
  ImageFiltering.jl whenever possible which grants a significant
  boost in performance.
* Surface correlation functions are much more accurate now. There is
  a new edge detection filter `Kernel5x5` which works for the most of
  cases and is default. An old kernel is now called `Kernel3x3`
  and is recommended for images taken with insufficient resolution
  (like `0.968 < lowfreq_energy_ratio(data) < 0.973`).
* `edgemode` parameter in surface function is renamed to `filter`
  and now must be of `EdgeFilter` type.

## Version 0.7.1

* Several bugs were fixed in the Map module.

## Version 0.7.0

* Improved tests and documentation.
* Edge detection mode must now be chosen with one of `Utilities.EdgeMode`
  types. These modes include a mode compatible with CorrelationTrackers.jl
  and a faster mode which works with the help of FFTW/CUFFT an is now
  default.
* `Map.mean_dir` was renamed to `Map.average_directions`.
* `Map.c2` is now about 2 times faster.
* Small performance improvements in `Directional.s2`.
* `Directional.c2` can now calculate the correlation function using
  FFT if a number of clusters is small enough.

## Version 0.6.1

* Compatibility release: Edge extraction on CPU behaves more like
  the previous method.

## Version 0.6.0

* Sobel operator is replaced with a more simple filter for edge
  detection.
* Significant speedup in Map.surfsurf and Map.surfvoid functions
  on GPU.
* lowfreq_energy_ratio and read_cuboid are moved to Utilities
  module.
* Map.l2 is currently removed for being obscure and inefficient.

## Version 0.5.0

* Tests consume less memory.
* Map module: add phase argument to all correlation functions to make
  API consistent with Directional module.
* Code cleanup.

## Version 0.4.5

* Improved functions from Map module.
