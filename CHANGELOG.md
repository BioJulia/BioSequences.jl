# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [3.0]
### Removed
- Removed `unsafe_setindex!`. Instead, use normal setindex with `@inbounds`.
- Removed minhashing functionality - see package MinHash.jl
- Removed composition functionality - see package Kmers.jl
- Removed ReferenceSequence functionality
- Removed demultiplexer functionality
- Removed kmer functionality - this is moved to Kmers.jl
- Removed VoidAlphabet and CharAlphabet
- Removed ConditionIterator

### Added
- Added type `LongSubSeq`, a view into a `LongSequence`.
- Added method `translate!(::LongAminoAcidSeq, ::LongNucleotideSeq; kwargs...)`
- Added method `join(::Type{T<:BioSeuence}, it)` to join an iterable of biosequences
to a new instance of T.
- Added method `join!(s::BioSequence, it)`, an in-place version of `join`

### Changed
- `LongSequence` is no longer copy-on-write. For views, use `LongSubSeq`.
- Renamed `LongAminoAcidSeq` -> `LongAA`, `LongNucleotideSeq` -> `LongNuc`
  `LongRNASeq` -> `LongRNA` and `LongDNASeq` -> `LongDNA`
- The interface for `Alphabet` and `BioSequence` is now more clearly defined, documented, and tested.
- The constructor `LongSequence{A}(::Integer)` has been removed in favor of `LongSequence{A}(undef, ::Integer)`.
- Biological sequences can no longer be converted to/from strings and vectors.
- Updated the element and substring search API to conform to `Base.find*` patterns.

## [2.0.1]
### Changed
- Fixed syntax errors where functions were marked with `@inbounds` instead of
  `@inline`.
  
## [2.0]
### Added
- New subtypes of Random.Sampler, SamplerUniform and SamplerWeighted.
- Random `LongSequence`s can now be created with `randseq`,
  optionally using a sampler to specify element distribution.
- All random `LongSequence` generator methods take an optional AbstractRNG
  argument.
- Add methods to `randseq` to optimize random generation of `NucleicAcid` or
  `AminoAcid` `LongSequence`s.
- BioGenerics is now a dependency - replaces BioCore.
- A `SkipmerFactory` iterator that allows iteration over the Skipmers in a 
  nucleotide sequence. A Skipmer is a `Mer` (see changed below), that is
  generated using a certain cyclic nucleotide sampling pattern.
  See [this paper](https://www.biorxiv.org/content/early/2017/09/19/179960.full.pdf+html)
  for more details.
- A `BigMer` parametric primitive type has been added, that has the same
  functionality as `Mer` (see changed section), but uses 128 bits instead of 64.
- An abstract parametric type called `AbstractMer` has been added to unify `Mer`
  and `BigMer`.
- Generators of bit-parallel iteration code have been introduced to help
  developers write bitparallel implementations of some methods. Counting GC
  content, matches and mismatches have been migrated to use these generators.
- Added `occursin` methods for exact matching.

### Changed
- The abstract `Sequence` type is now called `BioSequence{A}`.
- The type previously called `BioSequence{A}` is now `LongSequence{A}`.
- `Kmers` are now a parametric primitive type: `Mer{A<:NucleicAcidAlphabet{2},K}`.
- `unsafe_setindex!` has been made systematic for all `setindex` methods as a 
  way of bypassing all bound checking and `orphan!` calls.
- Kmer string literals have been updated, they are now `mer""` string literals,
  and they have a flag to enforce the type of `Mer` e.g.: `mer"ATCG"dna`,
  `mer"AUCG"rna`
- No longer use an old version of Twiddle and deprecated functions.
- Using `Base.count` with certain functions and sequence combinations dispatches
  to highly optimized bit-parallel implementations, falling back to a default
  naive counting loop by default for all other predicate-sequence combinations.
- No more implicit conversion from strings to biological sequences. The `Base.convert`
  methods have been renamed to `Base.parse` methods.

### Removed
- The FASTQ module.
- The FASTA module.
- The TwoBit module.
- The ABIF module.
- BioCore is no longer a dependency.
- Automa is no longer a dependency.

## [1.1.0]
### Changed
- Automatic conversion of `LongDNASeq` to `LongRNASeq` when translating
  sequences.
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
- A bug fix for `FASTA.Record` writing where the width parameter of a
  `FASTA.Writer` is less than or equal to zero.

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

[Unreleased]: https://github.com/BioJulia/BioSequences.jl/compare/v2.0.1...HEAD
[2.0.1]: https://github.com/BioJulia/BioSequences.jl/compare/v2.0.0...v2.0.1
[2.0.0]: https://github.com/BioJulia/BioSequences.jl/compare/v1.1.0...v2.0.0
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
