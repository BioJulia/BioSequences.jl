# Change Log

## [0.8.2](https://github.com/BioJulia/BioSequences.jl/tree/0.8.2) (2018-02-19)
[Full Changelog](https://github.com/BioJulia/BioSequences.jl/compare/v0.8.1...0.8.2)

**Fixed bugs:**

- Bugfix: Writing FASTA sequences when width \<= 0 [\#33](https://github.com/BioJulia/BioSequences.jl/pull/33) ([Ward9250](https://github.com/Ward9250))

**Closed issues:**

- Add split function [\#23](https://github.com/BioJulia/BioSequences.jl/issues/23)

**Merged pull requests:**

- Message printed when sequence is empty. [\#31](https://github.com/BioJulia/BioSequences.jl/pull/31) ([Ward9250](https://github.com/Ward9250))

## [v0.8.1](https://github.com/BioJulia/BioSequences.jl/tree/v0.8.1) (2017-11-10)
[Full Changelog](https://github.com/BioJulia/BioSequences.jl/compare/v0.8.0...v0.8.1)

**Implemented enhancements:**

- implement bit-parallel GC counter [\#29](https://github.com/BioJulia/BioSequences.jl/pull/29) ([bicycle1885](https://github.com/bicycle1885))

**Merged pull requests:**

- fix type definition keywords [\#28](https://github.com/BioJulia/BioSequences.jl/pull/28) ([bicycle1885](https://github.com/bicycle1885))
- Update documentation generation [\#26](https://github.com/BioJulia/BioSequences.jl/pull/26) ([Ward9250](https://github.com/Ward9250))

## [v0.8.0](https://github.com/BioJulia/BioSequences.jl/tree/v0.8.0) (2017-08-16)
[Full Changelog](https://github.com/BioJulia/BioSequences.jl/compare/v0.7.0...v0.8.0)

**Closed issues:**

- Drop Julia 0.5 support [\#14](https://github.com/BioJulia/BioSequences.jl/issues/14)

**Merged pull requests:**

- \[WIP\] Add unique sequence counting function [\#22](https://github.com/BioJulia/BioSequences.jl/pull/22) ([Ward9250](https://github.com/Ward9250))
- use Dict to store counts in Composition [\#21](https://github.com/BioJulia/BioSequences.jl/pull/21) ([bicycle1885](https://github.com/bicycle1885))
- add typemin and typemax for Kmer [\#20](https://github.com/BioJulia/BioSequences.jl/pull/20) ([bicycle1885](https://github.com/bicycle1885))
- Make minhash function for generic reader [\#19](https://github.com/BioJulia/BioSequences.jl/pull/19) ([kescobo](https://github.com/kescobo))
- update doctests for Julia 0.6 [\#18](https://github.com/BioJulia/BioSequences.jl/pull/18) ([bicycle1885](https://github.com/bicycle1885))
- add position-weighted matrix search [\#16](https://github.com/BioJulia/BioSequences.jl/pull/16) ([bicycle1885](https://github.com/bicycle1885))

## [v0.7.0](https://github.com/BioJulia/BioSequences.jl/tree/v0.7.0) (2017-07-28)
[Full Changelog](https://github.com/BioJulia/BioSequences.jl/compare/v0.6.3...v0.7.0)

**Merged pull requests:**

- \[Review\] Drop Julia 0.5 [\#17](https://github.com/BioJulia/BioSequences.jl/pull/17) ([Ward9250](https://github.com/Ward9250))

## [v0.6.3](https://github.com/BioJulia/BioSequences.jl/tree/v0.6.3) (2017-07-06)
[Full Changelog](https://github.com/BioJulia/BioSequences.jl/compare/v0.6.1...v0.6.3)

**Merged pull requests:**

- use IterTools instead of Iterators [\#13](https://github.com/BioJulia/BioSequences.jl/pull/13) ([bicycle1885](https://github.com/bicycle1885))
- Added diag\_val methods [\#12](https://github.com/BioJulia/BioSequences.jl/pull/12) ([Ward9250](https://github.com/Ward9250))

## [v0.6.1](https://github.com/BioJulia/BioSequences.jl/tree/v0.6.1) (2017-06-20)
[Full Changelog](https://github.com/BioJulia/BioSequences.jl/compare/v0.6.0...v0.6.1)

**Closed issues:**

- BioSequences.jl 0.6 [\#10](https://github.com/BioJulia/BioSequences.jl/issues/10)
- We need a codon type for work with coding sequences. [\#3](https://github.com/BioJulia/BioSequences.jl/issues/3)

**Merged pull requests:**

- Counting hotfix [\#11](https://github.com/BioJulia/BioSequences.jl/pull/11) ([Ward9250](https://github.com/Ward9250))

## [v0.6.0](https://github.com/BioJulia/BioSequences.jl/tree/v0.6.0) (2017-06-13)
[Full Changelog](https://github.com/BioJulia/BioSequences.jl/compare/v0.5.0...v0.6.0)

**Merged pull requests:**

- print info only when interactive [\#9](https://github.com/BioJulia/BioSequences.jl/pull/9) ([bicycle1885](https://github.com/bicycle1885))
- An eachkmer iterator fix [\#8](https://github.com/BioJulia/BioSequences.jl/pull/8) ([Ward9250](https://github.com/Ward9250))
- Added convenience method for removing gaps from a sequence. [\#7](https://github.com/BioJulia/BioSequences.jl/pull/7) ([Ward9250](https://github.com/Ward9250))
- use dirname\(@\_\_FILE\_\_\) instead of Pkg.dir in one more location [\#6](https://github.com/BioJulia/BioSequences.jl/pull/6) ([tkelman](https://github.com/tkelman))
- \[WIP\] Improve Julia 0.6 compatibility [\#5](https://github.com/BioJulia/BioSequences.jl/pull/5) ([bicycle1885](https://github.com/bicycle1885))
- check gap existence before making a kmer [\#4](https://github.com/BioJulia/BioSequences.jl/pull/4) ([bicycle1885](https://github.com/bicycle1885))

## [v0.5.0](https://github.com/BioJulia/BioSequences.jl/tree/v0.5.0) (2017-06-07)
[Full Changelog](https://github.com/BioJulia/BioSequences.jl/compare/v0.5.0-prerelease...v0.5.0)

**Closed issues:**

- Error tagging new release [\#1](https://github.com/BioJulia/BioSequences.jl/issues/1)

## [v0.5.0-prerelease](https://github.com/BioJulia/BioSequences.jl/tree/v0.5.0-prerelease) (2017-06-06)


\* *This Change Log was automatically generated by [github_changelog_generator](https://github.com/skywinder/Github-Changelog-Generator)*