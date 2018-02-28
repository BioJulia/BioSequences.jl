# BioSequences.jl

| **Release**                                                     | **Documentation**                                                               | **Maintainers**                             |
|:---------------------------------------------------------------:|:-------------------------------------------------------------------------------:|:-------------------------------------------:|
| [![][release-img]][release-url] [![][license-img]][license-url] | [![][docs-stable-img]][docs-stable-url] [![][docs-latest-img]][docs-latest-url] | ![][maintainer-a-img] ![][maintainer-b-img] |


## Description

BioSequences.jl provides DNA, RNA and amino acid sequence data types for the
julia language, with a comprehensive set of methods for common operations and
IO of major sequence data formats.   


## Installation

Install BioSequences from the Julia REPL:

```julia
julia> Pkg.add("BioSequences")
```

If you are interested in the cutting edge of the development, please check out
the master branch to try new features before release.


## Testing

BioSequences.jl is tested against Julia `0.6` and current `0.7-dev` on Linux, OS X, and Windows.

| **PackageEvaluator**                                            | **Latest Build Status**                                                                                |
|:---------------------------------------------------------------:|:------------------------------------------------------------------------------------------------------:|
| [![][pkg-0.6-img]][pkg-0.6-url] [![][pkg-0.7-img]][pkg-0.7-url] | [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] [![][codecov-img]][codecov-url]        |


## Contributing and Questions

We appreciate contributions from users including reporting bugs, fixing issues,
improving performance and adding new features.
Please go to the [contributing section of the documentation](biojulia.net/Contributing/latest)
for more information.

If you have a question about
contributing or using this package, you are encouraged to use the
[Bio category of the Julia discourse
site](https://discourse.julialang.org/c/domain/bio).

[release-img]: https://img.shields.io/github/release/BioJulia/BioSequences.jl.svg
[release-url]: https://github.com/BioJulia/BioSequences.jl/releases/latest

[license-img]: https://img.shields.io/badge/license-MIT-green.svg
[license-url]: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://biojulia.github.io/BioSequences.jl/latest
[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://biojulia.github.io/BioSequences.jl/stable

[maintainer-a-img]: https://img.shields.io/badge/BioJulia%20Maintainer-Ward9250-orange.svg
[maintainer-b-img]: https://img.shields.io/badge/BioJulia%20Maintainer-bicycle1885-orange.svg

[pkg-0.6-img]: https://pkg.julialang.org/badges/BioSequences_0.6.svg
[pkg-0.6-url]: https://pkg.julialang.org/detail/BioSequences
[pkg-0.7-img]: https://pkg.julialang.org/badges/BioSequences_0.7.svg
[pkg-0.7-url]: https://pkg.julialang.org/detail/BioSequences

[travis-img]: https://img.shields.io/travis/BioJulia/BioSequences.jl/master.svg?label=Linux+/+macOS
[travis-url]: https://travis-ci.org/BioJulia/BioSequences.jl

[appveyor-img]: https://img.shields.io/appveyor/ci/BioJulia/BioSequences.jl/master.svg?label=Windows
[appveyor-url]: https://ci.appveyor.com/project/Ward9250/biosequences-jl/branch/master

[codecov-img]: https://codecov.io/gh/BioJulia/BioSequences.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/BioJulia/BioSequences.jl
