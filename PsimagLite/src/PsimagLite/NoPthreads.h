// SPDX-FileCopyrightText: Copyright (c) 2009-2013, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup PsimagLite */
/*@{*/

/*! \file NoPthreads.h
 *
 *  If you want to say no to pthreads
 *
 */
#ifndef NO_PTHREADS_HEADER_H
#define NO_PTHREADS_HEADER_H

#include <cassert>
#include <iostream>

namespace PsimagLite {
template <typename PthreadFunctionHolderType> class NoPthreads {

public:

	NoPthreads(SizeType npthreads = 1, int = 0)
	{
		std::cerr << "NoPthreads is deprecated, please use NoPthreadsNg\n";
		assert(npthreads == 1);
	}

	void loopCreate(SizeType total, PthreadFunctionHolderType& pfh)
	{
		pfh.thread_function_(0, total, total, 0);
	}

	String name() const { return "nopthreads"; }

	SizeType threads() const { return 1; }

	SizeType mpiProcs() const { return 1; }

}; // NoPthreads
} // namespace PsimagLite
/*@}*/
#endif
