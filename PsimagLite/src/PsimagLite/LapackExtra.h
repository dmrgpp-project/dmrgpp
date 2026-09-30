// BEGIN LICENSE BLOCK
/*
Copyright (c) 2009 , UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 1.0.0]

*********************************************************
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED.

Please see full open source license included in file LICENSE.
*********************************************************

*/
// END LICENSE BLOCK

#ifndef LAPACK_H_
#define LAPACK_H_
#include <complex>

// See PsimagLite/LAPACK.h for the rationale: when MKL's own LAPACK header is
// available, reuse its (const-correct, width-correct) declarations instead
// of re-declaring the same extern "C" symbols with a possibly-conflicting
// signature. Falls back to the hand-written declarations below for any other
// LAPACK implementation (reference LAPACK, OpenBLAS, etc.)
#if defined(__has_include)
#if __has_include(<mkl_lapack.h>)
#define PSIMAGLITE_LAPACKEXTRA_MKL 1
#ifndef MKL_Complex8
#define MKL_Complex8 std::complex<float>
#endif
#ifndef MKL_Complex16
#define MKL_Complex16 std::complex<double>
#endif
#include <mkl_lapack.h>
#endif
#endif

#if defined(PSIMAGLITE_LAPACKEXTRA_MKL)
using IntegerForLapackType = MKL_INT;
#elif !defined(PSI_LAPACK_64)
using IntegerForLapackType = int;
#else
using IntegerForLapackType = long int;
#endif

#ifndef PSIMAGLITE_LAPACKEXTRA_MKL
extern "C" void zheev_(char*,
                       char*,
                       IntegerForLapackType*,
                       std::complex<double>*,
                       IntegerForLapackType*,
                       double*,
                       std::complex<double>*,
                       IntegerForLapackType*,
                       double*,
                       IntegerForLapackType*);
extern "C" void cheev_(char*,
                       char*,
                       IntegerForLapackType*,
                       std::complex<float>*,
                       IntegerForLapackType*,
                       float*,
                       std::complex<float>*,
                       IntegerForLapackType*,
                       float*,
                       IntegerForLapackType*);
extern "C" void dsyev_(char*,
                       char*,
                       IntegerForLapackType*,
                       double*,
                       IntegerForLapackType*,
                       double*,
                       double*,
                       IntegerForLapackType*,
                       IntegerForLapackType*);
extern "C" void ssyev_(char*,
                       char*,
                       IntegerForLapackType*,
                       float*,
                       IntegerForLapackType*,
                       float*,
                       float*,
                       IntegerForLapackType*,
                       IntegerForLapackType*);

extern "C" void zgeev_(char*,
                       char*,
                       IntegerForLapackType*,
                       std::complex<double>*,
                       IntegerForLapackType*,
                       std::complex<double>*,
                       std::complex<double>*,
                       IntegerForLapackType*,
                       std::complex<double>*,
                       IntegerForLapackType*,
                       std::complex<double>*,
                       IntegerForLapackType*,
                       double*,
                       IntegerForLapackType*);

extern "C" void dgeev_(char*,
                       char*,
                       IntegerForLapackType*,
                       double*,
                       IntegerForLapackType*,
                       double*,
                       double*,
                       double*,
                       IntegerForLapackType*,
                       double*,
                       IntegerForLapackType*,
                       double*,
                       IntegerForLapackType*,
                       IntegerForLapackType*);
#endif // !PSIMAGLITE_LAPACKEXTRA_MKL
#endif
