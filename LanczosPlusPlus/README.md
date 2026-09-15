# Quick Start
## Licensing

The full software license for LanczosPlusPlus++ version 1.0.0
can be found in
file [`LICENSE`](../LICENSE) in the repository root.
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

## DISCLAIMER

THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT OWNER, CONTRIBUTORS, UNITED STATES GOVERNMENT,
OR THE UNITED STATES DEPARTMENT OF ENERGY BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED
STATES DEPARTMENT OF ENERGY, NOR THE COPYRIGHT OWNER, NOR
ANY OF THEIR EMPLOYEES, REPRESENTS THAT THE USE OF ANY
INFORMATION, DATA, APPARATUS, PRODUCT, OR PROCESS
DISCLOSED WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.


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

