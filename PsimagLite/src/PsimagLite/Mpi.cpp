// SPDX-FileCopyrightText: Copyright (c) 2009-2013, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup PsimagLite */
/*@{*/

/*! \file Mpi.cpp
 *
 */

#ifdef DMRGPP_USE_MPI
#include "MpiYes.cpp"
#else
#include "MpiNo.cpp"
#endif
