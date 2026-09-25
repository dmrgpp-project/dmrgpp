// SPDX-FileCopyrightText: Copyright (c) 2009-2011, 2013, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file CorrectionVectorFunction.h
 *
 *  This is an implementation of PRB 60, 335, Eq. (24)
 *
 */
#ifndef CORRECTION_V_FUNCTION_H
#define CORRECTION_V_FUNCTION_H
#include "ConjugateGradient.h"
#include <PsimagLite/FreqEnum.h>

namespace Dmrg {
template <typename MatrixType, typename InfoType> class CorrectionVectorFunction {

	using FieldType  = typename MatrixType::value_type;
	using VectorType = typename PsimagLite::Vector<FieldType>::Type;
	using RealType   = typename PsimagLite::Real<FieldType>::Type;

	class InternalMatrix {

	public:

		using value_type = FieldType;
		InternalMatrix(const MatrixType& m, const InfoType& info, RealType E0)
		    : m_(m)
		    , info_(info)
		    , E0_(E0)
		{
			if (info_.omega().first != PsimagLite::FreqEnum::REAL)
				throw PsimagLite::RuntimeError("Matsubara only with KRYLOV\n");
		}

		SizeType rows() const { return m_.rows(); }

		void matrixVectorProduct(VectorType& x, const VectorType& y) const
		{
			RealType   eta          = info_.eta();
			RealType   omegaMinusE0 = info_.omega().second + E0_;
			VectorType xTmp(x.size(), 0);
			m_.matrixVectorProduct(xTmp, y); // xTmp = Hy
			VectorType x2(x.size(), 0);
			m_.matrixVectorProduct(x2, xTmp); // x2 = H^2 y
			const RealType f1 = (-2.0);
			// this needs fixing
			// preferred:
			// x <= x2 + f1*omegaMinusE0*xTmp + (omegaMinusE0*omegaMinusE0 + eta*eta)*y;
			// equivalent
			for (SizeType i = 0; i < x.size(); ++i)
				x[i] = x2[i] + f1 * omegaMinusE0 * xTmp[i]
				    + (omegaMinusE0 * omegaMinusE0 + eta * eta) * y[i];

			x /= (-eta);
		}

	private:

		const MatrixType& m_;
		const InfoType&   info_;
		RealType          E0_;
	};

	using ConjugateGradientType = ConjugateGradient<InternalMatrix>;

public:

	CorrectionVectorFunction(const MatrixType& m, const InfoType& info, RealType E0)
	    : im_(m, info, E0)
	    , cg_(info.cgSteps(), info.cgEps())
	{ }

	void getXi(VectorType& result, const VectorType& sv) const
	{
		VectorType x0(result.size(), 0.0);

		result = x0; // initial ansatz
		cg_(result, im_, sv);
	}

private:

	InternalMatrix        im_;
	ConjugateGradientType cg_;
}; // class CorrectionVectorFunction
} // namespace Dmrg

/*@}*/
#endif // CORRECTION_V_FUNCTION_H
