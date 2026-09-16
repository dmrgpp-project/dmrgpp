# Quick Start
## Licensing

LanczosPlusPlus first-party software is licensed under the GNU General Public
License, version 3 or, at your option, any later version
(`GPL-3.0-or-later`). See the root [`LICENSE`](../LICENSE) file. This software
comes with absolutely no warranty; see the GNU General Public License for
details. Institutional notices are in [`NOTICE`](../NOTICE), and terms for
bundled third-party material are in
[`THIRD_PARTY_NOTICES`](../THIRD_PARTY_NOTICES).

LanczosPlusPlus is a free and open source implementation of the
Lanczos algorithm for models of strongly correlated electrons.
You are welcomed to use it and publish data
obtained with Lanczos++. If you do, please cite this
work .

## Code Integrity

Hash of the latest commit is also posted at
https://g1257.github.com/hashes.html

Latest commit should always be signed.
https://g1257.github.com/keys.html
## How To Cite This Work

@article{re:alvarez09,
author="G. Alvarez",
title="The density matrix renormalization group for strongly correlated electron
systems: A generic implementation",
journal="Computer Physics Communications",
volume="180",
pages="1572",
year="2009"}

And also:

@article{
re:webDmrgPlusPlus,
Author = {G. Alvarez},
Title = {DMRG++ Website},
Publisher = {\url{https://g1257.github.com/dmrgPlusPlus}} }

Building and Running Lanczos++
TBW.

## Required Software

- (required) GNU C++
- (required) The LAPACK and BLAS libraries
This library is available for most platforms.
The configure.pl script will ask for the LDFLAGS variable
to pass to the compiler/linker. If the linux platform was
chosen the default/suggested LDFLAGS will include -llapack.
If the osx platform was chosen the default/suggested LDFLAGS will
include  -framework Accelerate.
For other platforms the appropriate linker flags must be given.
More information on LAPACK is here: http://netlib.org/lapack/
- (required) PsimagLite. This is here \url{https://github.com/g1257/PsimagLite/}.
You can do \verb=git clone https://github.com/g1257/PsimagLite.git= in a separate directory
outside of the DMRG++ distribution. \verb=configure.pl= will ask you where you put it.
- (optional) make or gmake (only needed to use the Makefile)
- (optional) perl (only needed to run the configure.pl script)

## Website

See https://g1257.github.com/LanczosPlusPlus/

