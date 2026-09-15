// SPDX-FileCopyrightText: Copyright (c) 2009-2017, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup PsimagLite */
/*@{*/

/*! \file Parallelizer.h
 *
 */
#ifndef PARALLELIZER_H
#define PARALLELIZER_H
#include "CodeSectionParams.h"
#include "Concurrency.h"
#include "Map.h"
#include "Vector.h"
#include <stdexcept>

#define ActualPthreadsName PthreadsNg
#define ActualPthreadsN() "PthreadsNg.h"
#define ActualNoPthreadsName NoPthreadsNg
#define ActualNoPthreadsN() "NoPthreadsNg.h"

#ifdef USE_PTHREADS

#include ActualPthreadsN()

#else

#include ActualNoPthreadsN()

#endif // USE_PTHREADS

namespace PsimagLite {
template <typename InstanceType, typename LoadBalancerType = LoadBalancerDefault>
class Parallelizer

#ifdef USE_PTHREADS

    : public ActualPthreadsName<InstanceType, LoadBalancerType> {
	using BaseType = ActualPthreadsName<InstanceType, LoadBalancerType>;

#else
    : public ActualNoPthreadsName<InstanceType, LoadBalancerType> {
	using BaseType = ActualNoPthreadsName<InstanceType, LoadBalancerType>;
#endif

public:

	Parallelizer(const CodeSectionParams& cs)
	    : BaseType(cs)
	{ }

	Parallelizer(String codeSectionParams)
	    : BaseType(codeSectionParamss_[codeSectionParams])
	{ }

	static bool exists(String name)
	{
		return (codeSectionParamss_.find(name) != codeSectionParamss_.end());
	}

	static void push(String name, const CodeSectionParams& cs)
	{
		codeSectionParamss_[name] = cs;
	}

	static void clear() { codeSectionParamss_.clear(); }

private:

	static Map<String, CodeSectionParams>::Type codeSectionParamss_;
};
} // namespace PsimagLite

/*@}*/
#endif
