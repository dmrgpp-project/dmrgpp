// SPDX-FileCopyrightText: Copyright (c) 2014, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup PsimagLite */
/*@{*/

/*! \file FloatingPoint.h
 *
 */
#ifndef PSI_FLOATING_POINT_H
#define PSI_FLOATING_POINT_H

#ifndef NDEBUG
#include <fenv.h>
#endif

namespace PsimagLite {

class FloatingPoint {

public:

	static void enableExcept()
	{
#ifndef NDEBUG
#ifndef __APPLE__
		feenableexcept(FE_DIVBYZERO | FE_INVALID);
#endif
#endif
	}
}; // class FloatingPoint

} // namespace PsimagLite

/*@}*/
#endif
