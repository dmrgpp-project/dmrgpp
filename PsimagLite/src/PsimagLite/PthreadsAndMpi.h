// SPDX-FileCopyrightText: Copyright (c) 2009-2013, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup PsimagLite */
/*@{*/

/*! \file PthreadsAndMpi .h
 *
 *  A C++ pthreads and MPI class that implements the Concurrency interface
 *
 */
#ifndef PTHREADS_AND_MPI_H
#define PTHREADS_AND_MPI_H

#include "NotMpi.h"
#include "Pthreads.h"
#include <iostream>
#include <pthread.h>

namespace PsimagLite {
template <typename PthreadFunctionHolderType>
class PthreadsAndMpi : public Pthreads<PthreadFunctionHolderType> {

	using BaseType = Pthreads<PthreadFunctionHolderType>;

public:

	PthreadsAndMpi(SizeType npthreads, MPI::CommType comm = MPI::COMM_WORLD)
	    : BaseType(npthreads)
	    , nthreads_(npthreads)
	    , comm_(comm)
	{ }

	String name() const { return "pthreadsandmpi"; }

	SizeType threads() const { return nthreads_; }

	SizeType mpiProcs() const { return MPI::commSize(comm_); }

private:

	SizeType        nthreads_;
	MPI::CommType   comm_;
	pthread_mutex_t mutex_;

}; // PthreadsAndMpi class

} // namespace PsimagLite

/*@}*/
#endif
