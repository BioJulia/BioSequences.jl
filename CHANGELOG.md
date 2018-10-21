# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.1.0]
### Changed
- Automatic conversion of `DNASequence` to `RNASequence` when translating sequences.
- Add `alternative_start` keyword argument to translate().
- Add abstract type for kmer iterators.
- :racehorse: Faster kmer iteration.
- Fixed indexing in ABIF records.

## [1.0.0] - 2018-08-23
### Added
- Issue and PR templates.
- Code of Conduct and Contributing files.
- A changelog file.
- Support for julia v0.7 and v1.0.

### Removed
- :exclamation: Support for julia v0.6.

## [0.8.3] - 2018-02-28
### Changed
- Fix the `sequence` method so as the sequence type can be specified, allowing
  type-stable efficient code generation.

## [0.8.2] - 2018-02-19
### Changed
- A bug fix for `FASTA.Record` writing where the width parameter of a `FASTA.Writer`
is less than or equal to zero.

## [0.8.1] - 2017-11-10
### Changed
- Update documentation generation.
- Fixes to type definition keywords.
- Bit-parallel GC counting.

## [0.8.0] - 2017-08-16
### Added
- Position weight matrix search functionality.
- A generalised composition method.
- `typemin` and `typemax` methods for `Kmer` types. 

### Changed
- `MinHash` function now generalised to `Reader` types.
- Updates to doc tests.

## [0.7.0] - 2017-07-28
### Added
- Support for julia v0.6 only. 

### Removed
- :exclamation: Dropped support for julia v0.5.

## [0.6.3] - 2017-07-06
### Changed
- Iterators.jl is not longer used as a dependency in favour of Itertools.jl.

## [0.6.1] - 2017-06-20
### Changed
- Bug-fix for site-counting algorithm.

## [0.6.0] - 2017-06-14
### Added
- :arrow_up: Compatibility with julia v0.6.
- The `ungap` and `ungap!` methods, that are shorthand for filtering gaps from
  biological sequences.

### Changed
- Bug fixes for Kmer iteration that were caused by gaps in 4-bit encoded sequences.

## [0.5.0] - 2017-06-07
### Added
- All files pertaining to the old Bio.Seq module.


[Unreleased]: https://github.com/BioJulia/BioSequences.jl/compare/v1.1.0...HEAD
[1.1.0]: https://github.com/BioJulia/BioSequences.jl/compare/v1.0.0...v1.1.0
[1.0.0]: https://github.com/BioJulia/BioSequences.jl/compare/v0.8.3...v1.0.0
[0.8.3]: https://github.com/BioJulia/BioSequences.jl/compare/v0.8.2...v0.8.3
[0.8.2]: https://github.com/BioJulia/BioSequences.jl/compare/v0.8.1...v0.8.2
[0.8.1]: https://github.com/BioJulia/BioSequences.jl/compare/v0.8.0...v0.8.1
[0.8.0]: https://github.com/BioJulia/BioSequences.jl/compare/v0.7.0...v0.8.0
[0.7.0]: https://github.com/BioJulia/BioSequences.jl/compare/v0.6.3...v0.7.0
[0.6.3]: https://github.com/BioJulia/BioSequences.jl/compare/v0.6.1...v0.6.3
[0.6.1]: https://github.com/BioJulia/BioSequences.jl/compare/v0.6.0...v0.6.1
[0.6.0]: https://github.com/BioJulia/BioSequences.jl/compare/v0.5.0...v0.6.0
[0.5.0]: https://github.com/BioJulia/BioSequences.jl/tree/v0.5.0
