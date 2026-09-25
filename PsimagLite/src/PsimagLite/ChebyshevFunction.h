
// SPDX-FileCopyrightText: Copyright (c) 2011, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]

/** \ingroup PsimagLite */
/*@{*/

/*! \file ChebyshevFunction.h
 *
 * The Chebyshev function
 *
 */

#ifndef CHEBYSHEV_FUNCTION_H
#define CHEBYSHEV_FUNCTION_H
#include "TypeToString.h"
#include <iostream>

namespace PsimagLite {
template <typename RealType> class ChebyshevFunction {

public:

	RealType operator()(int m, const RealType& x) const
	{
		if (m == 0)
			return 1;

		if (m == 1)
			return x;

		if (m & 1) {
			int p = (m - 1) / 2;
			return (2 * this->operator()(p, x) * this->operator()(p + 1, x) - x);
		}

		int      pp  = m / 2;
		RealType tmp = this->operator()(pp, x);
		return (2 * tmp * tmp - 1);
	}
}; // class ChebyshevFunction

} // namespace PsimagLite
/*@}*/
#endif // CHEBYSHEV_FUNCTION_H
