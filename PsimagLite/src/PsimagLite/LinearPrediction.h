// SPDX-FileCopyrightText: Copyright (c) 2009 , UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup PsimagLite */
/*@{*/

/*! \file LinearPrediction.h
 *
 *  Extrapolating a "time" series
 *  see extrapolation.tex for more details
 */

#ifndef LINEAR_PREDICTION_H
#define LINEAR_PREDICTION_H
#include "BLAS.h"
#include "LAPACK.h"
#include "Matrix.h"
#include <complex>

extern "C"
{
#include <gsl/gsl_poly.h>
}

namespace PsimagLite {

template <typename FieldType> class LinearPrediction {
	using MatrixType = Matrix<FieldType>;

public:

	LinearPrediction(const typename Vector<FieldType>::Type& y, SizeType p)
	    : y_(y)
	    , p_(p)
	{
		SizeType ysize = y.size();
		if (ysize & 1)
			throw RuntimeError("LinearPrediction::ctor(...): data set must "
			                   "contain an even number of points\n");
		SizeType                         n = ysize / 2;
		MatrixType                       A(p, p);
		typename Vector<FieldType>::Type B(p);
		computeA(A, n);
		computeB(B, n);
		computeD(A, B);
	}

	const FieldType& operator()(SizeType i) const { return y_[i]; }

	void linearPredictionfunction(const typename Vector<FieldType>::Type& y, SizeType p)
	{
		SizeType ysize = y.size();
		if (ysize & 1)
			throw RuntimeError("LinearPrediction::ctor(...): data set must "
			                   "contain an even number of points\n");
		SizeType                         n = ysize / 2;
		MatrixType                       A(p, p);
		typename Vector<FieldType>::Type B(p);
		computeA(A, n);
		computeB(B, n);
		computeD(A, B);
	}

	void predict(SizeType p)
	{
		SizeType n = y_.size();
		for (SizeType i = n; i < n + p; i++) {
			FieldType sum = 0;
			for (SizeType j = 0; j < d_.size(); j++) {
				sum += d_[j] * y_[i - j - 1];
			}
			y_.push_back(sum);
		}
	}

	void predict2(SizeType p)
	{
		// Fix roots
		using Complex                           = std::complex<double>;
		SizeType                         twicep = 2 * p;
		SizeType                         pp1    = p + 1;
		typename Vector<FieldType>::Type nd(p), aa(pp1), zz(twicep);
		std::vector<Complex>             roots(p);
		std::vector<Complex>             ab(pp1);

		aa[p] = 1;
		for (SizeType j = 0; j < p; j++) {
			aa[j] = -d_[p - 1 - j];
		}

		gsl_poly_complex_workspace* w = gsl_poly_complex_workspace_alloc(pp1);
		gsl_poly_complex_solve(&aa[0], pp1, w, &zz[0]);
		gsl_poly_complex_workspace_free(w);

		for (SizeType j = 0; j < p; j++) {
			roots[j] = Complex(zz[2 * j], zz[2 * j + 1]);
		}

		for (SizeType j = 0; j < p; j++) {
			// Look for a root outside the unit circle, and put it
			// to 0
			if (abs(roots[j]) > 1.0) {
				roots[j] = 0.0; // roots[j]/abs(roots[j]);
			}
		}
		// Now reconstruct the polynomial coefficients
		ab[0] = -roots[0];
		ab[1] = 1.0;
		for (SizeType j = 1; j < p; j++) {
			ab[j + 1] = 1.0;
			for (SizeType i = j; i >= 1; i--) {
				ab[i] = ab[i - 1] - roots[j] * ab[i];
			}
			ab[0] = -roots[j] * ab[0];
		}
		for (SizeType j = 0; j < p; j++) {
			nd[p - 1 - j] = -PsimagLite::real(ab[j]);
		}

		SizeType n = y_.size();
		y_.resize(n + p);
		for (SizeType i = n; i < n + p; i++) {
			FieldType sum = 0;
			for (SizeType j = 0; j < d_.size(); j++) {
				sum += nd[j] * y_[i - j - 1];
			}
			y_[i] = sum;
		}
	}

private:

	//! Note: A and B cannot be const. here due to the ultimate
	//! call to BLAS::GEMV
	void computeD(MatrixType& A, typename Vector<FieldType>::Type& B)
	{
		SizeType                   p = B.size();
		typename Vector<int>::Type ipiv(p); // use signed integers here!!
		int                        info = 0;
		psimag::LAPACK::GETRF(p, p, &(A(0, 0)), p, &(ipiv[0]), info);

		typename Vector<FieldType>::Type work(2);
		int                              lwork = -1; // query mode
		psimag::LAPACK::GETRI(p, &(A(0, 0)), p, &(ipiv[0]), &(work[0]), lwork, info);
		lwork = static_cast<int>(work[0]);
		if (lwork <= 0)
			throw RuntimeError("LinearPrediction:: internal error\n");
		work.resize(lwork);
		// actual work:
		psimag::LAPACK::GETRI(p, &(A(0, 0)), p, &(ipiv[0]), &(work[0]), lwork, info);

		d_.resize(p);
		psimag::BLAS::GEMV('N', p, p, 1.0, &(A(0, 0)), p, &(B[0]), 1, 0.0, &(d_[0]), 1);
	}

	void computeA(MatrixType& A, SizeType n) const
	{
		SizeType p = A.rows();
		for (SizeType l = 0; l < p; l++) {
			for (SizeType j = 0; j < p; j++) {
				A(l, j) = 0;
				for (SizeType i = n; i < 2 * n; i++)
					A(l, j) += y_[i - l - 1] * y_[i - j - 1];
			}
		}
	}

	void computeB(typename Vector<FieldType>::Type& B, SizeType n) const
	{
		SizeType p = B.size();
		for (SizeType l = 0; l < p; l++) {
			B[l] = 0;
			for (SizeType i = n; i < 2 * n; i++)
				B[l] += y_[i - l - 1] * y_[i];
		}
	}

	typename Vector<FieldType>::Type y_;
	SizeType                         p_;
	typename Vector<FieldType>::Type d_;
}; // class LinearPrediction
} // namespace PsimagLite

/*@}*/
#endif // LINEAR_PREDICTION_H
