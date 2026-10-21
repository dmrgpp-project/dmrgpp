# DMRG++ Licensing and Third-Party Provenance Audit

**Audit date:** 2026-09-14

**Target license:** `GPL-3.0-or-later`

**Baseline:** commit `c3760bc5f8575c6c7b8982ad0d15e2e0918fe614`

**Scope:** All 1,383 paths in the baseline commit on branch `license_change`; untracked files, subsequently added audit records, build trees, fetched dependency caches, and local generated output were excluded.

## Result

The tracked-source audit is complete. It found 304 tracked files with at least one former-license fingerprint: 285 C/C++ source files and 19 legal, documentation, manual-input, or script files. No tracked file currently contains an SPDX declaration or a project GPL declaration.

No classified item found by this audit still blocks bulk relicensing or GPL-covered distribution. The three provenance issues found in the baseline were resolved after the audit by removing one third-party work and independently replacing two implementations, as recorded below.

## Reproducible inventory method

The audit uses only paths tracked by the recorded baseline commit. The former-license manifest is the union of case-insensitive matches for:

- `UT Battelle Open Source Software License`
- `All rights reserved`
- `DE-AC05-00OR22725`
- `THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS`
- `DISCLOSED WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS`
- `full open source license included in file LICENSE`

Reproduce the core counts with:

```sh
git ls-tree -r --name-only -z c3760bc5f8575c6c7b8982ad0d15e2e0918fe614 | tr '\0' '\n' | wc -l
git grep -Il -i \
  -e 'UT Battelle Open Source Software License' \
  -e 'All rights reserved' \
  -e 'DE-AC05-00OR22725' \
  -e 'THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS' \
  -e 'DISCLOSED WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS' \
  -e 'full open source license included in file LICENSE' \
  c3760bc5f8575c6c7b8982ad0d15e2e0918fe614
git grep -Il -E 'SPDX-(License-Identifier|FileCopyrightText)' \
  c3760bc5f8575c6c7b8982ad0d15e2e0918fe614
git grep -Il -i -E 'GNU (General Public License|GPL)|GPL-[0-9]|GPLv?[0-9]' \
  c3760bc5f8575c6c7b8982ad0d15e2e0918fe614
git grep -In -i -E \
  'adapted from|translated from|copied from|parts taken from|sample code|original source|non-commercial|permission is granted|permission to use|inspired by|source:' \
  c3760bc5f8575c6c7b8982ad0d15e2e0918fe614
```

Observed component counts:

| Component | Former-license files |
|---|---:|
| `(root)` | 1 |
| `LanczosPlusPlus` | 32 |
| `PsimagLite` | 95 |
| `dmrg` | 174 |
| `doc` | 1 |
| `scripts` | 1 |

The narrower `All rights reserved` search matches 301 files. The complete union adds `LanczosPlusPlus/README.md`, `PsimagLite/README.md`, and `dmrg/KronUtil/README`.

## Legal and documentation surfaces

Tracked legal files:

- `LICENSE`
- `dmrg/LICENSE`
- `LanczosPlusPlus/LICENSE`
- `PsimagLite/LICENSE`

`PsimagLite/LICENSE:69-86` uniquely centralizes the Andrei Alexandrescu Loki notice, while `PsimagLite/loki/Sequence.h:2-10` carries Peter Kümmel's separate notice. Both are preserved in `THIRD_PARTY_NOTICES` before any nested-license removal.

Stale human-facing declarations occur at:

- `README.md:16-20`
- `cincuenta/README.md:2-11`
- `LanczosPlusPlus/README.md:2-7,42-65`
- `PsimagLite/README.md:6-7,16-39`
- `dmrg/KronUtil/README:6-33`

Manual-generation inputs copy or embed the former license:

- `LanczosPlusPlus/doc/manual.ptex:87-119,132,196-199`
- `PsimagLite/doc/manual.ptex:92-123,136,243-246`
- `doc/manual.ptex:87-119,978-981`; the final directive references nonexistent `../src/LICENSE` relative to `doc/`.

No tracked test fixture or golden output contains a former-license fingerprint.

## Runtime surfaces

Compiled former-license banners are defined at:

- `dmrg/ProgramGlobals.cpp:84-101`
- `LanczosPlusPlus/src/lanczos.cpp:3-17`
- `LanczosPlusPlus/src/ed.cpp:4-18`

Observed print paths are:

- `dmrg/dmrg.cpp:14-23,115-121`
- `dmrg/toolboxdmrg.cpp:133-138`
- `LanczosPlusPlus/src/lanczos.cpp:175-179`
- `LanczosPlusPlus/src/ed.cpp:102-106` (not currently a CMake executable target)

`dmrg/introspect.cpp:15-24` defines an unused banner-printing helper. `dmrg/manyOmegas.cpp:76-82` prints provenance but not the former license banner.

## Third-party path-level classification

### Compatible; terms retained

**Loki Library — permissive terms**

Andrei Alexandrescu notice:

- `PsimagLite/loki/NullType.h:2-13`
- `PsimagLite/loki/TypeManip.h:2-13`
- `PsimagLite/loki/TypeTraits.h:2-13`
- `PsimagLite/loki/Typelist.h:2-13`
- `PsimagLite/loki/TypelistMacros.h:2-13`

Peter Kümmel notice:

- `PsimagLite/loki/Sequence.h:2-10`

The permission texts are GPL-compatible and are retained verbatim in `THIRD_PARTY_NOTICES`.

**René Nyffenegger Base64 — permissive zlib-style terms**

- `PsimagLite/src/PsimagLite/PsiBase64.h:3-33`
- `PsimagLite/src/PsimagLite/PsiBase64.cpp:3-34`
- `PsimagLite/examples/base64test.cpp:3-34`

All three carry the permission text and modification statement. Their terms are retained verbatim in `THIRD_PARTY_NOTICES`.

**Stack Overflow SFINAE example — CC BY-SA 4.0**

- `PsimagLite/examples/sfinaeHasMember.cpp:4-22`
- Source: <https://stackoverflow.com/a/257382>, Nicola Bonelli (user 19630)
- The Stack Exchange API returned `content_license: CC BY-SA 4.0` on 2026-09-14.

The repository version is adapted. Creative Commons lists GPLv3 as a one-way BY-SA-compatible license for adaptations of CC BY-SA 4.0 material, subject to attribution, change-notice, and other adapter requirements: <https://creativecommons.org/share-your-work/cclicenses/compatible-licenses/>. `THIRD_PARTY_NOTICES` records the creator, source, license, and changes. Preserve the source comment or replace it with an equivalent in-file pointer when headers are updated.

### Unresolved or incompatible; blocks migration

None of the items identified in the baseline audit remains unresolved.

### Resolved after the baseline audit

**Stack Overflow permutation parity implementation — independently replaced**

- Baseline implementation: `dmrg/Engine/Observer.h:601-625`
- Baseline source: <https://stackoverflow.com/a/20718842>, mikera (user 214010)
- Baseline license: CC BY-SA 3.0

This finding was resolved on 2026-09-14 by removing the adapted `swapCountSmall` implementation and its source citation. `dmrg/Engine/PermutationParity.h` now contains an independently authored inversion-parity implementation, which `Observer.h` uses. Regression coverage in `dmrg/Engine/tests/test_PermutationParity.cpp` verifies basic cycles, identity and reversed permutations, sizes beyond machine-word bit masks, randomized transposition sequences, and inverse permutations. Because the third-party implementation is no longer distributed, its section was removed from `THIRD_PARTY_NOTICES`.

**Farid Khoury sparse matrix multiplication — independently replaced**

- Baseline implementation: `PsimagLite/src/PsimagLite/CrsMatrix.h:980-1027`
- Baseline citation: <http://web.maths.unsw.edu.au/~farid/Papers/Hons/node23.html>
- Archived Fortran implementation: <https://web.archive.org/web/20160304000000id_/http://web.maths.unsw.edu.au/~farid/Papers/Hons/node41.html>

This finding was resolved on 2026-09-14 by removing the translated implementation and its source citation and replacing it with an independently authored row-local associative-accumulator implementation. Regression coverage in `PsimagLite/src/tests/test_CrsMatrixMultiply.cpp` verifies rectangular and random products, canonical column ordering, exact-zero cancellation, empty rows, identity multiplication, mixed scalar types, dimension validation, and output aliasing. Because the third-party implementation is no longer distributed, its section was removed from `THIRD_PARTY_NOTICES`.

**Wikimedia under-construction icon — removed**

- Baseline paths: `doc/Under_construction_icon-blue.pdf` and `doc/Under_construction_icon-blue.eps`
- Baseline embedding/credit: `doc/manual.ptex:130-133`
- Source: <https://commons.wikimedia.org/wiki/File:Under_construction_icon-blue.svg>, author Dsmurat

The Wikimedia source page states LGPL version 2.1 or later and requires attribution. This finding was resolved on 2026-09-14 by deleting both tracked derivatives and removing the corresponding inclusion and credit from `doc/manual.ptex`. Because the third-party work is no longer distributed, its section was removed from `THIRD_PARTY_NOTICES`.

### Attribution/reference only; no separate license boundary found

- `PsimagLite/src/PsimagLite/Geometry/LongRange.h:323-330` credits an algorithmic idea without asserting copied source.
- `cincuenta/TestSuite/gbek_reference/README.md:3-22` describes an independent implementation of published equations.
- `cincuenta/TestSuite/gbek_reference/fetch_arxiv_figures.sh:2-9` downloads paper figures only into an untracked, ignored directory; distributing generated reports containing them requires a separate rights review.
- `PsimagLite/src/PsimagLite/MersenneTwister.h:33-47`, `PsimagLite/src/PsimagLite/RungeKutta.h:78`, `scripts/honeycomb.pl:2-3`, `scripts/honeycombObservables.pl:2-3`, `PsimagLite/src/PsimagLite/GemmR.h:8`, and `dmrg/KronUtil/README:2-7` contain first-party contributor credits covered by the approved first-party scope; preserve the names.
- `cincuenta/src/ImpuritySolverNeqGBEK.h:31-58` cites a paper and equations without asserting copied sample code.
- `cincuenta/src/NeqBathDecomposition.h:365-366` identifies a port of this repository's own tracked Python reference implementation.

The prior untracked notice entries for `cincuenta/src/KadanoffBaym.h` and `dmrg/Engine/ClebschGordan.h` were removed because neither path exists or is tracked in the audited tree.

## Dependency and distribution boundary

Tracked CMake installs targets but no legal files:

- `dmrg/CMakeLists.txt:62-65`
- `LanczosPlusPlus/src/CMakeLists.txt:24`
- `cincuenta/src/CMakeLists.txt:10`
- `PsimagLite/src/CMakeLists.txt:112-114`
- `dmrg/KronUtil/CMakeLists.txt:4-6`

No tracked `install(FILES ...)`, CPack rule, or package rule installs `LICENSE`, `NOTICE`, or `THIRD_PARTY_NOTICES`.

Catch2, Kokkos, and Kokkos Kernels are fetched dependencies; MPI, HDF5, Boost, Threads, BLAS, LAPACK, GSL, and optional MAGMA are system dependencies. They are not vendored tracked source and are outside this source-tree notice inventory. Binary distributors must still account for the terms of dependencies they redistribute.

## Complete former-license manifest (304 tracked files)

<details>
<summary>Show paths</summary>

- `LICENSE`
- `LanczosPlusPlus/LICENSE`
- `LanczosPlusPlus/README.md`
- `LanczosPlusPlus/doc/manual.ptex`
- `LanczosPlusPlus/src/Engine/BasisBase.h`
- `LanczosPlusPlus/src/Engine/DefaultSymmetry.h`
- `LanczosPlusPlus/src/Engine/Engine.h`
- `LanczosPlusPlus/src/Engine/InputCheck.h`
- `LanczosPlusPlus/src/Engine/InternalProductOnTheFly.h`
- `LanczosPlusPlus/src/Engine/InternalProductStored.h`
- `LanczosPlusPlus/src/Engine/LanczosGlobals.h`
- `LanczosPlusPlus/src/Engine/LanczosModelBase.hpp`
- `LanczosPlusPlus/src/Engine/Partitions.h`
- `LanczosPlusPlus/src/Engine/ReducedDensityMatrix.h`
- `LanczosPlusPlus/src/Engine/ReflectionSymmetry.h`
- `LanczosPlusPlus/src/Engine/TranslationSymmetry.h`
- `LanczosPlusPlus/src/Models/FeBasedSc/BasisFeAsBasedSc.h`
- `LanczosPlusPlus/src/Models/FeBasedSc/BasisFeAsSpinOrbit.h`
- `LanczosPlusPlus/src/Models/FeBasedSc/BasisOneSpinFeAs.h`
- `LanczosPlusPlus/src/Models/FeBasedSc/FeBasedSc.h`
- `LanczosPlusPlus/src/Models/FeBasedSc/ParametersModelFeAs.h`
- `LanczosPlusPlus/src/Models/FermionSpinless/ParametersFermionSpinless.hh`
- `LanczosPlusPlus/src/Models/Heisenberg/ParametersHeisenberg.h`
- `LanczosPlusPlus/src/Models/HubbardOneOrbital/ParametersModelHubbard.h`
- `LanczosPlusPlus/src/Models/Immm/BasisImmm.h`
- `LanczosPlusPlus/src/Models/Immm/BasisOneSpinImmm.h`
- `LanczosPlusPlus/src/Models/Immm/Immm.h`
- `LanczosPlusPlus/src/Models/Immm/ParametersImmm.h`
- `LanczosPlusPlus/src/Models/Kitaev/ParametersKitaev.h`
- `LanczosPlusPlus/src/Models/TjMultiOrb/ParametersTjMultiOrb.h`
- `LanczosPlusPlus/src/configure.pl`
- `LanczosPlusPlus/src/ed.cpp`
- `LanczosPlusPlus/src/lanczos.cpp`
- `PsimagLite/LICENSE`
- `PsimagLite/PsimagDoc/README`
- `PsimagLite/README.md`
- `PsimagLite/TestSuite/README`
- `PsimagLite/doc/manual.ptex`
- `PsimagLite/examples/affinityTest.cpp`
- `PsimagLite/examples/akimaSpline.cpp`
- `PsimagLite/examples/combineContinuedFraction.cpp`
- `PsimagLite/examples/continuedFraction.cpp`
- `PsimagLite/examples/continuedFractionCollection.cpp`
- `PsimagLite/examples/kernelPolynomial.cpp`
- `PsimagLite/examples/linearPrediction.cpp`
- `PsimagLite/examples/loadImbalance.cpp`
- `PsimagLite/examples/nested.cpp`
- `PsimagLite/examples/newconfigure.pl`
- `PsimagLite/examples/sparseSolverTest.cpp`
- `PsimagLite/examples/threads.cpp`
- `PsimagLite/lib/configure.pl`
- `PsimagLite/lib/newconfigure.pl`
- `PsimagLite/src/PsimagLite/AST/AdditionalFunctions.h`
- `PsimagLite/src/PsimagLite/AST/ExpressionForAST.h`
- `PsimagLite/src/PsimagLite/AST/Node.h`
- `PsimagLite/src/PsimagLite/AST/PlusMinusMultiplyDivide.h`
- `PsimagLite/src/PsimagLite/AST/Tree.h`
- `PsimagLite/src/PsimagLite/Ainur/configure.pl`
- `PsimagLite/src/PsimagLite/AllocatorCpu.h`
- `PsimagLite/src/PsimagLite/AlmostEqual.h`
- `PsimagLite/src/PsimagLite/ApplicationInfo.h`
- `PsimagLite/src/PsimagLite/ChebyshevFunction.h`
- `PsimagLite/src/PsimagLite/ChebyshevFunctionExplicit.h`
- `PsimagLite/src/PsimagLite/ChebyshevSerializer.h`
- `PsimagLite/src/PsimagLite/ChebyshevSolver.h`
- `PsimagLite/src/PsimagLite/Complex.h`
- `PsimagLite/src/PsimagLite/Concurrency.h`
- `PsimagLite/src/PsimagLite/ContinuedFraction.h`
- `PsimagLite/src/PsimagLite/ContinuedFractionCollection.h`
- `PsimagLite/src/PsimagLite/CrsMatrix.h`
- `PsimagLite/src/PsimagLite/Fermi.h`
- `PsimagLite/src/PsimagLite/FloatingPoint.h`
- `PsimagLite/src/PsimagLite/Geometry/Geometry.h`
- `PsimagLite/src/PsimagLite/Geometry/GeometryBase.h`
- `PsimagLite/src/PsimagLite/Geometry/GeometryDca.h`
- `PsimagLite/src/PsimagLite/Geometry/GeometryDirection.h`
- `PsimagLite/src/PsimagLite/Geometry/GeometryTerm.h`
- `PsimagLite/src/PsimagLite/Geometry/Honeycomb.h`
- `PsimagLite/src/PsimagLite/Geometry/KTwoNiFFour.h`
- `PsimagLite/src/PsimagLite/Geometry/Ladder.h`
- `PsimagLite/src/PsimagLite/Geometry/LadderBath.h`
- `PsimagLite/src/PsimagLite/Geometry/LadderX.h`
- `PsimagLite/src/PsimagLite/Geometry/LongChain.h`
- `PsimagLite/src/PsimagLite/Geometry/LongRange.h`
- `PsimagLite/src/PsimagLite/Geometry/Star.h`
- `PsimagLite/src/PsimagLite/GslWrapper.h`
- `PsimagLite/src/PsimagLite/InputCheckBase.h`
- `PsimagLite/src/PsimagLite/InputNg.h`
- `PsimagLite/src/PsimagLite/Io/IoNg.h`
- `PsimagLite/src/PsimagLite/Io/IoSimple.h`
- `PsimagLite/src/PsimagLite/LanczosCore.h`
- `PsimagLite/src/PsimagLite/LanczosVectors.h`
- `PsimagLite/src/PsimagLite/LapackExtra.h`
- `PsimagLite/src/PsimagLite/LineMarker.h`
- `PsimagLite/src/PsimagLite/LinearPrediction.h`
- `PsimagLite/src/PsimagLite/Map.h`
- `PsimagLite/src/PsimagLite/Matrix.h`
- `PsimagLite/src/PsimagLite/MatrixSolverBase.hh`
- `PsimagLite/src/PsimagLite/MemoryCpu.h`
- `PsimagLite/src/PsimagLite/MemoryUsage.h`
- `PsimagLite/src/PsimagLite/Mpi.cpp`
- `PsimagLite/src/PsimagLite/MpiNo.h`
- `PsimagLite/src/PsimagLite/MpiYes.cpp`
- `PsimagLite/src/PsimagLite/MpiYes.h`
- `PsimagLite/src/PsimagLite/NoPthreads.h`
- `PsimagLite/src/PsimagLite/NoPthreadsNg.h`
- `PsimagLite/src/PsimagLite/NotMpi.h`
- `PsimagLite/src/PsimagLite/Options.h`
- `PsimagLite/src/PsimagLite/PackIndices.h`
- `PsimagLite/src/PsimagLite/Parallelizer.h`
- `PsimagLite/src/PsimagLite/ParametersForSolver.h`
- `PsimagLite/src/PsimagLite/Permutations.h`
- `PsimagLite/src/PsimagLite/Profiling.h`
- `PsimagLite/src/PsimagLite/ProgressIndicator.cpp`
- `PsimagLite/src/PsimagLite/ProgressIndicator.h`
- `PsimagLite/src/PsimagLite/Pthreads.h`
- `PsimagLite/src/PsimagLite/PthreadsAndMpi.h`
- `PsimagLite/src/PsimagLite/PthreadsNg.h`
- `PsimagLite/src/PsimagLite/Random48.h`
- `PsimagLite/src/PsimagLite/RandomForTests.h`
- `PsimagLite/src/PsimagLite/RungeKutta.h`
- `PsimagLite/src/PsimagLite/SampleCRSMatrix.h`
- `PsimagLite/src/PsimagLite/SparseRow.h`
- `PsimagLite/src/PsimagLite/SparseRowCached.h`
- `PsimagLite/src/PsimagLite/SparseVector.h`
- `PsimagLite/src/PsimagLite/Stack.h`
- `PsimagLite/src/PsimagLite/TridiagonalMatrix.h`
- `PsimagLite/src/PsimagLite/Vector.h`
- `dmrg/Engine/ApplyOperatorExpression.h`
- `dmrg/Engine/ApplyOperatorLocal.h`
- `dmrg/Engine/Basis.h`
- `dmrg/Engine/BasisWithOperators.h`
- `dmrg/Engine/BlockDiagonalMatrix.h`
- `dmrg/Engine/Checkpoint.h`
- `dmrg/Engine/Cloner.h`
- `dmrg/Engine/ConjugateGradient.h`
- `dmrg/Engine/CorrectionVectorFunction.h`
- `dmrg/Engine/CorrectionVectorSkeleton.h`
- `dmrg/Engine/CorrelationsSkeleton.h`
- `dmrg/Engine/DensityMatrixBase.h`
- `dmrg/Engine/DensityMatrixLocal.h`
- `dmrg/Engine/DensityMatrixSvd.h`
- `dmrg/Engine/Diagonalization.h`
- `dmrg/Engine/DiskStackNg.h`
- `dmrg/Engine/DmrgSerializer.h`
- `dmrg/Engine/DmrgSolver.h`
- `dmrg/Engine/FermionSign.h`
- `dmrg/Engine/FourPointCorrelations.h`
- `dmrg/Engine/HamiltonianConnection.h`
- `dmrg/Engine/HamiltonianSymmetryLocal.h`
- `dmrg/Engine/InputCheck.h`
- `dmrg/Engine/JmPairs.h`
- `dmrg/Engine/LeftRightSuper.h`
- `dmrg/Engine/LineChangerLinear.h`
- `dmrg/Engine/Link.h`
- `dmrg/Engine/MatrixOrIdentity.h`
- `dmrg/Engine/MatrixVectorBase.h`
- `dmrg/Engine/MatrixVectorKron/ArrayOfMatStruct.h`
- `dmrg/Engine/MatrixVectorKron/GenIjPatch.h`
- `dmrg/Engine/MatrixVectorKron/InitKronBase.h`
- `dmrg/Engine/MatrixVectorKron/InitKronHamiltonian.h`
- `dmrg/Engine/MatrixVectorKron/KronConnections.h`
- `dmrg/Engine/MatrixVectorKron/KronMatrix.h`
- `dmrg/Engine/MatrixVectorKron/MatrixVectorKron.h`
- `dmrg/Engine/MatrixVectorOnTheFly.h`
- `dmrg/Engine/MatrixVectorStored.h`
- `dmrg/Engine/MettsCollapse.h`
- `dmrg/Engine/MettsParams.h`
- `dmrg/Engine/MettsSerializer.h`
- `dmrg/Engine/MettsStochastics.h`
- `dmrg/Engine/ModelBase.h`
- `dmrg/Engine/ModelCommon.h`
- `dmrg/Engine/ModelHelperLocal.h`
- `dmrg/Engine/ModelSelector.h`
- `dmrg/Engine/MultiPointCorrelations.h`
- `dmrg/Engine/ObservableLibrary.h`
- `dmrg/Engine/ObservablesOnePointInSitu.h`
- `dmrg/Engine/Observer.h`
- `dmrg/Engine/ObserverHelper.h`
- `dmrg/Engine/OnePointCorrelations.h`
- `dmrg/Engine/Operator.h`
- `dmrg/Engine/Operators.h`
- `dmrg/Engine/Parallel2PointCorrelations.h`
- `dmrg/Engine/Parallel4PointDs.h`
- `dmrg/Engine/ParallelDensityMatrix.h`
- `dmrg/Engine/ParallelTriDiag.h`
- `dmrg/Engine/ParametersDmrgSolver.h`
- `dmrg/Engine/ProgramGlobals.h`
- `dmrg/Engine/Recovery.h`
- `dmrg/Engine/Runner.h`
- `dmrg/Engine/SpinSquared.h`
- `dmrg/Engine/SpinSquaredHelper.h`
- `dmrg/Engine/Su2Related.h`
- `dmrg/Engine/TargetHelper.h`
- `dmrg/Engine/TargetParamsBase.h`
- `dmrg/Engine/TargetParamsCommon.h`
- `dmrg/Engine/TargetParamsCorrection.h`
- `dmrg/Engine/TargetParamsCorrectionVector.h`
- `dmrg/Engine/TargetParamsDynamic.h`
- `dmrg/Engine/TargetParamsGroundState.h`
- `dmrg/Engine/TargetParamsTimeStep.h`
- `dmrg/Engine/TargetParamsTimeVectors.h`
- `dmrg/Engine/TargetQuantumElectrons.h`
- `dmrg/Engine/TargetingBase.h`
- `dmrg/Engine/TargetingCVEvolution.h`
- `dmrg/Engine/TargetingChebyshev.h`
- `dmrg/Engine/TargetingCommon.h`
- `dmrg/Engine/TargetingCorrection.h`
- `dmrg/Engine/TargetingCorrectionVector.h`
- `dmrg/Engine/TargetingDynamic.h`
- `dmrg/Engine/TargetingExpression.h`
- `dmrg/Engine/TargetingGroundState.h`
- `dmrg/Engine/TargetingMetts.h`
- `dmrg/Engine/TargetingRixsDynamic.h`
- `dmrg/Engine/TargetingRixsStatic.h`
- `dmrg/Engine/TargetingTimeStep.h`
- `dmrg/Engine/TimeSerializer.h`
- `dmrg/Engine/TimeVectorsBase.h`
- `dmrg/Engine/TimeVectorsChebyshev.h`
- `dmrg/Engine/TimeVectorsKrylov.h`
- `dmrg/Engine/TimeVectorsRungeKutta.h`
- `dmrg/Engine/TimeVectorsSuzukiTrotter.h`
- `dmrg/Engine/ToolBox.h`
- `dmrg/Engine/TridiagRixsStatic.h`
- `dmrg/Engine/Truncation.h`
- `dmrg/Engine/TwoPointCorrelations.h`
- `dmrg/Engine/Utils.h`
- `dmrg/Engine/VectorWithOffsets.h`
- `dmrg/Engine/VerySparseMatrix.h`
- `dmrg/Engine/Wft/ParallelWftMany.h`
- `dmrg/Engine/Wft/ParallelWftOne.h`
- `dmrg/Engine/Wft/ParallelWftSu2.h`
- `dmrg/Engine/Wft/WaveFunctionTransfBase.h`
- `dmrg/Engine/Wft/WaveFunctionTransfFactory.h`
- `dmrg/Engine/Wft/WaveFunctionTransfLocal.h`
- `dmrg/Engine/Wft/WaveStructSvd.h`
- `dmrg/KronUtil/README`
- `dmrg/LICENSE`
- `dmrg/Models/ExtendedHubbard1Orb/ExtendedHubbard1Orb.h`
- `dmrg/Models/FeAsBasedScExtended/FeAsBasedScExtended.h`
- `dmrg/Models/FeAsModel/HilbertSpaceFeAs.h`
- `dmrg/Models/FeAsModel/ModelFeBasedSc.h`
- `dmrg/Models/FeAsModel/ParametersModelFeAs.h`
- `dmrg/Models/FermionSpinless/FermionSpinless.h`
- `dmrg/Models/FermionSpinless/HilbertSpaceFermionSpinless.h`
- `dmrg/Models/FermionSpinless/ParametersFermionSpinless.h`
- `dmrg/Models/GaugeSpin/GaugeSpin.h`
- `dmrg/Models/GaugeSpin/ParametersGaugeSpin.h`
- `dmrg/Models/Graphene/Graphene.h`
- `dmrg/Models/Graphene/ParametersGraphene.h`
- `dmrg/Models/Heisenberg/HeisenbergMix.h`
- `dmrg/Models/Heisenberg/ModelHeisenberg.h`
- `dmrg/Models/Heisenberg/ParametersModelHeisenberg.h`
- `dmrg/Models/HeisenbergAncillaC/HeisenbergAncillaC.h`
- `dmrg/Models/HeisenbergAncillaC/ParametersHeisenbergAncillaC.h`
- `dmrg/Models/HolsteinSpinlessThin/HolsteinSpinlessThin.h`
- `dmrg/Models/HolsteinThin/HolsteinThin.h`
- `dmrg/Models/HubbardAncilla/HubbardAncilla.h`
- `dmrg/Models/HubbardAncilla/ParametersHubbardAncilla.h`
- `dmrg/Models/HubbardAncillaExtended/HubbardAncillaExtended.h`
- `dmrg/Models/HubbardAncillaExtended/ParametersHubbardAncillaExtended.h`
- `dmrg/Models/HubbardHolstein/HilbertSpaceHubbardHolstein.h`
- `dmrg/Models/HubbardHolstein/HubbardHolstein.h`
- `dmrg/Models/HubbardHolstein/ParametersHubbardHolstein.h`
- `dmrg/Models/HubbardHolsteinSpinless/HilbertSpaceHubbardHolsteinSpinless.h`
- `dmrg/Models/HubbardHolsteinSpinless/HubbardHolsteinSpinless.h`
- `dmrg/Models/HubbardHolsteinSpinless/ParametersHubbardHolsteinSpinless.h`
- `dmrg/Models/HubbardMultiBand/ModelHubbardMultiBand.h`
- `dmrg/Models/HubbardMultiBand/ParametersHubbardMultiBand.h`
- `dmrg/Models/HubbardOneBand/HilbertSpaceHubbard.h`
- `dmrg/Models/HubbardOneBand/ModelHubbard.h`
- `dmrg/Models/HubbardOneBand/ParametersModelHubbard.h`
- `dmrg/Models/Immm/HilbertSpaceImmm.h`
- `dmrg/Models/Immm/Immm.h`
- `dmrg/Models/Immm/ParametersImmm.h`
- `dmrg/Models/IsingMultiOrb/HilbertSpaceIsingMultiOrb.h`
- `dmrg/Models/IsingMultiOrb/ModelIsingMultiOrb.h`
- `dmrg/Models/IsingMultiOrb/ParametersModelIsingMultiOrb.h`
- `dmrg/Models/Kitaev/Kitaev.h`
- `dmrg/Models/Kitaev/ParametersKitaev.h`
- `dmrg/Models/LiouvillianHeisenberg/LiouvillianHeisenberg.hh`
- `dmrg/Models/LiouvillianHeisenberg/ParamsLiouvillianHeisenberg.hh`
- `dmrg/Models/SpinOrbital/ParametersSpinOrbital.h`
- `dmrg/Models/SpinOrbital/SpinOrbitalModel.h`
- `dmrg/Models/Su3/ParametersSu3.h`
- `dmrg/Models/Su3/Su3Model.h`
- `dmrg/Models/SuperExtendedHubbard1Orb/SuperExtendedHubbard1Orb.h`
- `dmrg/Models/SuperHubbardExtended/SuperHubbardExtended.h`
- `dmrg/Models/TjAncillaC/ParametersTjAncillaC.h`
- `dmrg/Models/TjAncillaC/TjAncillaC.h`
- `dmrg/Models/TjAncillaC2/TjAncillaC2.h`
- `dmrg/Models/TjAncillaG/TjAncillaG.h`
- `dmrg/Models/TjAnisotropic/ParametersModelTjAnisotropic.h`
- `dmrg/Models/TjAnisotropic/TjAnisotropic.h`
- `dmrg/Models/TjMultiOrb/ParametersModelTjMultiOrb.h`
- `dmrg/Models/TjMultiOrb/TjMultiOrb.h`
- `dmrg/Models/UlsOsu/ParametersModelUlsOsu.h`
- `dmrg/Models/UlsOsu/UlsOsu.h`
- `dmrg/ProgramGlobals.cpp`
- `dmrg/Su2Related.cpp`
- `dmrg/Utils.cpp`
- `dmrg/run_drivers.pl`
- `doc/manual.ptex`
- `scripts/input.pl`

</details>

## Task 2 disposition

The tracked-source inventory and classification are complete. All identified blockers have been resolved by removal or independent replacement. `THIRD_PARTY_NOTICES` has been reconciled to the current tracked tree and records the compatible third-party terms that remain distributed.
