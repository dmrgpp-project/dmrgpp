// SPDX-FileCopyrightText: Copyright (c) 2009-2013, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]

// A class to profile a scope of code
#ifndef PROFILING_H_
#define PROFILING_H_

#include "MemoryUsage.h"
#include "ProgressIndicator.h"
#include "PsimagLite.h"

#include <Kokkos_Profiling_ScopedRegion.hpp>

#include <iostream>

namespace PsimagLite {

class Profiling {

	double diff(double x1, double x2) const { return x1 - x2; }

public:

	Profiling(String caller, String additional, std::ostream& os)
	    : progressIndicator_(caller)
	    , memoryUsage_("/proc/self/stat")
	    , isDead_(false)
	    , os_(os)
	    , kokkos_region_(caller + (additional.empty() ? "" : (' ' + additional)))
	{
		OstringStream msg(std::cout.precision());
		msg() << "starting clock " << additional;
		progressIndicator_.printline(msg, os);
	}

	Profiling(String caller, std::ostream& os)
	    : Profiling(caller, {}, os)
	{ }

	Profiling(const Profiling&) = delete;

	Profiling& operator=(const Profiling&) = delete;

	~Profiling() { killIt(""); }

	void end(String message) { killIt(message); }

private:

	void killIt(String message)
	{
		if (isDead_)
			return;
		OstringStream msg(std::cout.precision());
		msg() << "stopping clock " << message;
		progressIndicator_.printline(msg, os_);
		isDead_ = true;
	}

	ProgressIndicator               progressIndicator_;
	MemoryUsage                     memoryUsage_;
	bool                            isDead_;
	std::ostream&                   os_;
	Kokkos::Profiling::ScopedRegion kokkos_region_;
}; // Profiling
} // namespace PsimagLite

#endif
