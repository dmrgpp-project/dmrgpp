
// SPDX-FileCopyrightText: Copyright (c) 2009 , UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]

/** \ingroup DMRG */
/*@{*/

/*! \file LineMarker.h
 *
 * Writes a line to disk, intended as a mark
 */

#ifndef LINE_MARKER_H
#define LINE_MARKER_H
#include "AllocatorCpu.h"
#include <iostream>

namespace PsimagLite {

class LineMarker {
public:

	LineMarker(const String& name)
	    : name_(name + "=0")
	{ }

	template <typename IoInputType>
	LineMarker(IoInputType& io, const String& name, SizeType level = 0)
	{
		name_      = name + "=0";
		SizeType x = 0; // bogus
		io.readline(x, name_, level);
	}

	template <typename IoOutputType> void save(IoOutputType& io) const { io.printline(name_); }

private:

	String name_;

}; // class LineMarker
} // namespace PsimagLite
/*@}*/
#endif // LINE_MARKER_H
