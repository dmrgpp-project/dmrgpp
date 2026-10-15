// SPDX-FileCopyrightText: Copyright (c) 2009, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 2.0.0]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file Runner.h
 *
 *  A class to help repeating a task multiple times
 *
 */
#ifndef RUNNER_HEADER_H
#define RUNNER_HEADER_H

#include "String.h"
#include <PsimagLite/TypeToString.h>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

namespace Dmrg {
class Runner {
public:

	Runner(const PsimagLite::String& exec,
	       const PsimagLite::String& inputRoot,
	       const PsimagLite::String& ext)
	    : exec_(exec)
	    , inputRoot_(inputRoot)
	    , ext_(ext)
	{ }

	void operator()(SizeType i) const
	{
		PsimagLite::String s = exec_ + " " + inputRoot_ + ttos(i) + ext_;
		system(s.c_str());
	}

private:

	PsimagLite::String exec_, inputRoot_, ext_;
}; // class Runner
} // namespace Dmrg

/*@}*/
#endif // RUNNER_HEADER_H
