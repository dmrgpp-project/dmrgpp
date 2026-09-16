// SPDX-FileCopyrightText: Copyright (c) 2009 , UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]

#ifndef LAPACK_H_
#define LAPACK_H_
#include <complex>

#ifndef PSI_LAPACK_64
using IntegerForLapackType = int;
#else
using IntegerForLapackType = long int;
#endif

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
                       double*,
                       IntegerForLapackType*);
#endif
