BioSequences.jl v0.8.0 Release Notes
====================================


Feature additions
-----------------

* Position weight matrix search ([#16]).
* Generalised composition method ([#22]).
* `typemin` and `typemax` methods have been added for `Kmer` types ([#20]).


Patches and fixes
-----------------

* Dictionaries now used to store counts in `Composition` types ([#21]).
* MinHash function is now generalised to `Reader` types ([#19]).
* Doctests have been updated to julia v0.6 ([#18]).

[#22]: https://github.com/BioJulia/BioSequences.jl/pull/22
[#21]: https://github.com/BioJulia/BioSequences.jl/pull/21
[#20]: https://github.com/BioJulia/BioSequences.jl/pull/20
[#19]: https://github.com/BioJulia/BioSequences.jl/pull/19
[#18]: https://github.com/BioJulia/BioSequences.jl/pull/18
[#16]: https://github.com/BioJulia/BioSequences.jl/pull/16
