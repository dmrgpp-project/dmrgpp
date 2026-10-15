
# Quick Start

## Introduction

PsimagLite first-party software is licensed under the GNU General Public
License, version 3 or, at your option, any later version
(`GPL-3.0-or-later`). See the root [`LICENSE`](../LICENSE) file. This software
comes with absolutely no warranty; see the GNU General Public License for
details. Institutional notices are in [`NOTICE`](../NOTICE), and terms for
bundled third-party material are in
[`THIRD_PARTY_NOTICES`](../THIRD_PARTY_NOTICES).

Please cite PsimagLite if you base any scientific
publication on this software. Citation should read:
G. Alvarez, (2011), PsimagLite (version 1.0)
[computer software], Oak Ridge National Laboratory.

## Description

PsimagLite is a collection of C++ classes that are common to
 codes for the simulation of strongly correlated electrons.
PsimagLite is inspired in T.S.'s Psimag software (but PsimagLite is not a fork of Psimag).

The reason for PsimagLite is to share code among different applications.
Applications that depend on PsimagLite are:
SpinPhononFermion, DMRG++, Lanczos++, FreeFermions, GpusDoneRight, BetheAnsatz

## Code integrity

Hash of the latest commit is also posted at
https://g1257.github.com/hashes.html

Latest commit should always be signed.
Keys at https://g1257.github.com/keys.html

## Building
git clone https://github.com/g1257/PsimagLite
cd PsimagLite/
git checkout features
git pull origin features
mkdir build
cd build
cmake ..
make -j number_here

