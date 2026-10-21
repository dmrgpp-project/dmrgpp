// SPDX-FileCopyrightText: Copyright (c) 2009-2013, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 2.]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup PsimagLite */
/*@{*/

/*! \file MpiNo.h
 *
 * Do NOT include this file directly, use Mpi.h instead
 *
 */
#ifndef MPINO_HEADER_H
#define MPINO_HEADER_H
#include "Vector.h"
#include "loki/TypeTraits.h"
#include <algorithm>
#include <stdexcept>

namespace PsimagLite {

namespace MPI {

	using CommType = int;
	extern int COMM_WORLD;
	extern int SUM;

	void init(int*, char***);

	void finalize();

	bool hasMpi();

	void info(std::ostream&);

	void version(std::ostream&);

	SizeType commSize(CommType);

	SizeType commRank(CommType);

	int barrier(CommType);

	template <typename T> void bcast(T&, int = 0, CommType = COMM_WORLD) { }

	template <typename T> void recv(T&, int, int, CommType = COMM_WORLD) { }

	template <typename T> void send(T&, int, int, CommType = COMM_WORLD) { }

	template <typename T> void pointByPointGather(T&, int = 0, CommType = COMM_WORLD) { }

	template <typename T> void reduce(T&, int = 0, int = 0, int = 0) { }

	template <typename T> void allReduce(T&) { }

} // namespace MPI

} // namespace PsimagLite

/*@}*/
#endif
