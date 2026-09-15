// SPDX-FileCopyrightText: Copyright (c) 2009-2013, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]

/** \ingroup PsimagLite */
/*@{*/

/*!
 *
 *
 */
#ifndef MEMORY_CPU_H
#define MEMORY_CPU_H

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace PsimagLite {

class MemoryCpu {

public:

	void deallocate(void* p)
	{
		assert(p);
		free(p);
		std::ostringstream msg;
		msg << "// Freed " << p;
		std::cout << msg.str() << "\n";
		p = 0;
	}

	void* allocate(size_t x)
	{
		void*              p = malloc(x);
		std::ostringstream msg;
		msg << "// Allocated starting at " << p << " " << x << " bytes";
		std::cout << msg.str() << "\n";
		return p;
	}

}; // class MemoryCpu

MemoryCpu globalMemoryCpu;
} // namespace PsimagLite

/*@}*/
#endif // MEMORY_CPU_H
