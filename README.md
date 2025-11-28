# <img src="./sticker.svg" width="30%" align="right" /> BioSequences

[![Latest Release](https://img.shields.io/github/release/BioJulia/BioSequences.jl.svg)](https://github.com/BioJulia/BioSequences.jl/releases/latest)
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://biojulia.github.io/BioSequences.jl/stable)
[![Pkg Status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

## Description

BioSequences provides data types and methods for common operations with
biological sequences, including DNA, RNA, and amino acid sequences.


## Installation

You can install BioSequences from the julia
REPL. Press `]` to enter pkg mode, and enter the following:

```julia 
add BioSequences
```
## Example
```julia
using BioSequences

# Create a DNA sequence
seq = DNASequence("ACGTACGT")
println(seq)

# Get the reverse complement
revcomp = reverse_complement(seq)
println(revcomp)
```

If you are interested in the cutting edge of the development, please check out
the master branch to try new features before release.


## Testing

BioSequences is tested against Julia `1.X` on Linux, OS X, and Windows.

[![](https://codecov.io/gh/BioJulia/BioSequences.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/BioJulia/BioSequences.jl)
[![Downstream](https://github.com/BioJulia/BioSequences.jl/actions/workflows/Downstream.yml/badge.svg)](https://github.com/BioJulia/BioSequences.jl/actions/workflows/Downstream.yml)

Run tests from the Julia REPL:

```julia
using Pkg
Pkg.test("BioSequences")
```
## Contributing

We appreciate contributions from users including reporting bugs, fixing
issues, improving performance and adding new features.

Take a look at the [contributing files](https://github.com/BioJulia/Contributing) for detailed contributor and maintainer guidelines, and code of conduct.

Steps to contribute:

1. Fork the repository
2. Clone your fork locally
3. Make your changes
4. Commit & push to your fork
5. Open a Pull Request


## Backers & Sponsors

Thank you to all our backers and sponsors!

## Questions?

If you have a question about contributing or using BioJulia software, come
on over and chat to us on [the Julia Slack workspace](https://julialang.org/slack/), or you can try the
[Bio category of the Julia discourse site](https://discourse.julialang.org/c/domain/bio).
