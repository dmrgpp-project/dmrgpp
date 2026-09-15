// SPDX-FileCopyrightText: Copyright (c) 2009 , UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup PsimagLite */
/*@{*/

/*! \file AlmostEqual.h
 *
 * Almost equal
 */

#ifndef ALMOST_EQUAL_H
#define ALMOST_EQUAL_H
#include "Matrix.h"

namespace PsimagLite {
template <typename RealType> class AlmostEqual {
public:

	AlmostEqual(const RealType& eps)
	    : eps_(eps)
	{ }

	//		template<typename SomeType>
	//		bool operator()(
	//			const Matrix<SomeType>& a,
	//			const Matrix<SomeType>& b)
	//		{
	//			for (SizeType i=0;i<a.n_row();i++) {
	//                for (SizeType j=0;j<a.n_col();j++) {
	//                    RealType x = PsimagLite::norm(a(i,j)-b(i,j));
	//                    if (x>eps_) {
	//						std::cerr<<"a("<<i<<","<<j<<")="<<a(i,j);
	//						std::cerr<<"
	// b("<<i<<","<<j<<")="<<b(i,j)<<"\n";
	// throw RuntimeError("almostEqual\n");
	//					}
	//				}
	//			}
	//			return true;
	//		}

private:

	const RealType& eps_;
}; // class AlmostEqual
} // namespace PsimagLite

/*@}*/
#endif // ALMOST_EQUAL_H
