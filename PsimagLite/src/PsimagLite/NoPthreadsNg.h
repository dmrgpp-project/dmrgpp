// SPDX-FileCopyrightText: Copyright (c) 2009-2017, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup PsimagLite */
/*@{*/

/*! \file NoPthreadsNg .h
 *
 *  A C++ NoPthreadsNg class that implements the Concurrency interface
 *
 */
#ifndef PSI_NOPTHREADS_NG_H
#define PSI_NOPTHREADS_NG_H
#include "CodeSectionParams.h"
#include "LoadBalancerDefault.h"

namespace PsimagLite {

template <typename PthreadFunctionHolderType, typename LoadBalancerType = LoadBalancerDefault>
class NoPthreadsNg {

public:

	using VectorSizeType = LoadBalancerDefault::VectorSizeType;

	NoPthreadsNg(const CodeSectionParams& cs)
	{
		if (cs.npthreads != 1)
			throw PsimagLite::RuntimeError("NoPthreadsNg: ctor with threads != 1\n");
	}

	bool affinities() const { return false; }

	size_t stackSize() const { return 0; }

	// no weights, no balancer ==> create weights, set all weigths to 1,
	// delegate
	void loopCreate(PthreadFunctionHolderType& pfh)
	{
		LoadBalancerType* loadBalancer = new LoadBalancerType(pfh.tasks(), 1);
		loopCreate(pfh, *loadBalancer);
		delete loadBalancer;
		loadBalancer = 0;
	}

	// weights, no balancer ==> create balancer with weights ==> delegate
	void loopCreate(PthreadFunctionHolderType& pfh, const VectorSizeType& weights)
	{
		LoadBalancerType* loadBalancer = new LoadBalancerType(weights, 1);
		loopCreate(pfh, *loadBalancer);
		delete loadBalancer;
		loadBalancer = 0;
	}

	// balancer (includes weights)
	void loopCreate(PthreadFunctionHolderType& pfh, const LoadBalancerType&)
	{
		SizeType ntasks = pfh.tasks();
		for (SizeType taskNumber = 0; taskNumber < ntasks; ++taskNumber) {
			pfh.doTask(taskNumber, 0);
		}
	}

	String name() const { return "NoPthreadsNg"; }

	SizeType threads() const { return 1; }

	SizeType mpiProcs() const { return 1; }

	void setAffinities(bool) { }
}; // NoPthreadsNg class

} // namespace PsimagLite

/*@}*/
#endif
