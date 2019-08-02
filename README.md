# <img src="./sticker.svg" width="30%" align="right" /> BioSequences

[![Latest Release](https://img.shields.io/github/release/BioJulia/BioSequences.jl.svg)](https://github.com/BioJulia/BioSequences.jl/releases/latest)
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE) 
[![Stable documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://biojulia.github.io/BioSequences.jl/stable)
[![Latest documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://biojulia.github.io/BioSequences.jl/latest/)
![Lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)
[![Chat](https://img.shields.io/gitter/room/BioJulia/BioSequences.svg)](https://gitter.im/BioJulia/BioSequences.jl)


## Description

BioSequences provides data types and methods for common operations with 
biological sequences, including DNA, RNA, and amino acid sequences.
It also provides I/O for common sequence file formats including FASTA, FASTQ,
2bit and more. 


## Installation

Releases of BioSequences version 2.0 and above are registered and made available
to install through BioJulia's package registry. Julia's package manager only
monitors the "General" package repository by default. So before you start, you
should tell julia about the existence of the BioJulia package registry.

Start a julia terminal, hit the `]` key to enter pkg mode (you should see the
prompt change from `julia>` to `pkg>` ), then enter the following command:

```julia
registry add https://github.com/BioJulia/BioJuliaRegistry.git
```

After you've added the registry, you can install BioSequences from the julia
REPL. Press `]` to enter pkg mode again, and enter the following:

```julia
add BioSequences
```

If you are interested in the cutting edge of the development, please check out
the master branch to try new features before release.


## Testing

BioSequences is tested against Julia `1.X` on Linux, OS X, and Windows.

[![](https://travis-ci.org/BioJulia/BioSequences.jl.svg?branch=master)](https://travis-ci.org/BioJulia/BioSequences.jl)
[![](https://ci.appveyor.com/api/projects/status/1vdxlfv7yk9c1kfb/branch/master?svg=true)](https://ci.appveyor.com/project/BenJWard/biosequences-jl/branch/master)
[![](https://codecov.io/gh/BioJulia/BioSequences.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/BioJulia/BioSequences.jl)|


## Contributing

We appreciate contributions from users including reporting bugs, fixing
issues, improving performance and adding new features.

Take a look at the [contributing files](https://github.com/BioJulia/Contributing)
detailed contributor and maintainer guidelines, and code of conduct.



### Financial contributions

We also welcome financial contributions in full transparency on our
[open collective](https://opencollective.com/biojulia).
Anyone can file an expense. If the expense makes sense for the development
of the community, it will be "merged" in the ledger of our open collective by
the core contributors and the person who filed the expense will be reimbursed.


## Backers & Sponsors

Thank you to all our backers and sponsors!

Love our work and community? [Become a backer](https://opencollective.com/biojulia#backer).

[![backers](https://opencollective.com/biojulia/backers.svg?width=890)](https://opencollective.com/biojulia#backers)

Does your company use BioJulia? Help keep BioJulia feature rich and healthy by
[sponsoring the project](https://opencollective.com/biojulia#sponsor)
Your logo will show up here with a link to your website.

[![](https://opencollective.com/biojulia/sponsor/0/avatar.svg)](https://opencollective.com/biojulia/sponsor/0/website)
[![](https://opencollective.com/biojulia/sponsor/1/avatar.svg)](https://opencollective.com/biojulia/sponsor/1/website)
[![](https://opencollective.com/biojulia/sponsor/2/avatar.svg)](https://opencollective.com/biojulia/sponsor/2/website)
[![](https://opencollective.com/biojulia/sponsor/3/avatar.svg)](https://opencollective.com/biojulia/sponsor/3/website)
[![](https://opencollective.com/biojulia/sponsor/4/avatar.svg)](https://opencollective.com/biojulia/sponsor/4/website)
[![](https://opencollective.com/biojulia/sponsor/5/avatar.svg)](https://opencollective.com/biojulia/sponsor/5/website)
[![](https://opencollective.com/biojulia/sponsor/6/avatar.svg)](https://opencollective.com/biojulia/sponsor/6/website)
[![](https://opencollective.com/biojulia/sponsor/7/avatar.svg)](https://opencollective.com/biojulia/sponsor/7/website)
[![](https://opencollective.com/biojulia/sponsor/8/avatar.svg)](https://opencollective.com/biojulia/sponsor/8/website)
[![](https://opencollective.com/biojulia/sponsor/9/avatar.svg)](https://opencollective.com/biojulia/sponsor/9/website)


## Questions?

If you have a question about contributing or using BioJulia software, come
on over and chat to us on [Gitter](https://gitter.im/BioJulia/General), or you can try the
[Bio category of the Julia discourse site](https://discourse.julialang.org/c/domain/bio).