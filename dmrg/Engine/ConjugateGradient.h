// SPDX-FileCopyrightText: Copyright (c) 2009-2014, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [by G.A., Oak Ridge National Laboratory]


/** \ingroup DMRG */
/*@{*/

/*! \file ConjugateGradient.h
 *
 *  impl. of the conjugate gradient method
 *
 */
#ifndef CONJ_GRAD_H
#define CONJ_GRAD_H

#include <PsimagLite/Matrix.h>
#include <PsimagLite/ProgressIndicator.h>
#include <PsimagLite/Vector.h>

namespace Dmrg {

template <typename MatrixType> class ConjugateGradient {
	using FieldType  = typename MatrixType::value_type;
	using VectorType = typename PsimagLite::Vector<FieldType>::Type;
	using RealType   = typename PsimagLite::Real<FieldType>::Type;

public:

	ConjugateGradient(SizeType max, RealType eps)
	    : progress_("ConjugateGradient")
	    , max_(max)
	    , eps_(eps)
	{ }

	//! A and b, the result x, and also the initial solution x0
	void operator()(VectorType& x, const MatrixType& A, const VectorType& b) const
	{
		VectorType v = multiply(A, x);
		VectorType p(b.size());
		VectorType rprev(b.size());
		VectorType rnext;
		for (SizeType i = 0; i < rprev.size(); i++) {
			rprev[i] = b[i] - v[i];
			p[i]     = rprev[i];
		}

		SizeType k = 0;
		while (k < max_) {
			VectorType    tmp         = multiply(A, p);
			FieldType     scalarrprev = scalarProduct(rprev, rprev);
			FieldType     val         = scalarrprev / scalarProduct(p, tmp);
			v <= x + val* p;
			x = v;
			v <= rprev - val* tmp;
			rnext = v;
			if (PsimagLite::norm(rnext) < eps_)
				break;
			val = scalarProduct(rnext, rnext) / scalarrprev;
			v <= rnext - val* p;
			p     = v;
			rprev = rnext;
			k++;
		}

		PsimagLite::OstringStream                     msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg << "Finished after " << k << " steps out of " << max_;
		msg << " requested eps= " << eps_;
		RealType finalEps = PsimagLite::norm(rnext);
		msg << " actual eps= " << finalEps;
		progress_.printline(msgg, std::cout);

		if (finalEps <= eps_)
			return;

		PsimagLite::OstringStream                     msgg2(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg2 = msgg2();
		msg2 << "WARNING: actual eps " << finalEps
		     << " greater than requested eps= " << eps_;
		progress_.printline(msgg2, std::cout);
	}

private:

	FieldType scalarProduct(const VectorType& v1, const VectorType& v2) const
	{
		FieldType sum = 0;
		for (SizeType i = 0; i < v1.size(); i++)
			sum += PsimagLite::conj(v1[i]) * v2[i];
		return sum;
	}

	VectorType multiply(const MatrixType& A, const VectorType& v) const
	{
		VectorType y(A.rows(), 0);
		A.matrixVectorProduct(y, v);
		return y;
	}

	PsimagLite::ProgressIndicator progress_;
	SizeType                      max_;
	RealType                      eps_;
}; // class ConjugateGradient

} // namespace Dmrg

/*@}*/
#endif // CONJ_GRAD_H
