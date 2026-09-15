// SPDX-FileCopyrightText: Copyright (c) 2009 , UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]
// [by G.A., Oak Ridge National Laboratory]


/** \ingroup PsimagLite */
/*@{*/

/*! \file PackIndices.h
 *
 * PAcking and unpacking of indices
 */

#ifndef PACK_INDICES_H
#define PACK_INDICES_H
#include <cstdlib>

namespace PsimagLite {
class PackIndices {

public:

	PackIndices(SizeType n)
	    : n_(n)
	{ }

	//! given ind and n, get x and y such that ind = x + y*n
	void unpack(SizeType& x, SizeType& y, SizeType ind) const
	{
		// y  = ind/n;
		// x = ind - y*n;
		// x= ind % n;
		ldiv_t q = std::ldiv(ind, n_);
		y        = q.quot;
		x        = q.rem;
	}

	SizeType
	pack(SizeType x0, SizeType x1, const Vector<SizeType>::Type& permutationInverse) const
	{
		assert(x0 + n_ * x1 < permutationInverse.size());
		return permutationInverse[x0 + n_ * x1];
	}

private:

	SizeType n_;
}; // class PackIndices
} // namespace PsimagLite

/*@}*/
#endif // PACK_INDICES_H
