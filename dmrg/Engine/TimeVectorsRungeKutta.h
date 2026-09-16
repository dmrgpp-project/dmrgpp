// SPDX-FileCopyrightText: Copyright (c) 2009-2012, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 2.0.0]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file TimeVectorsRungeKutta.h
 *
 *
 */

#ifndef TIME_VECTORS_RUNGE_KUTTA
#define TIME_VECTORS_RUNGE_KUTTA
#include "TimeVectorsBase.h"
#include <PsimagLite/RungeKutta.h>
#include <iostream>

namespace Dmrg {

template <typename RealType> RealType minusOneOrMinusI(RealType) { return -1; }

template <typename RealType> std::complex<RealType> minusOneOrMinusI(std::complex<RealType>)
{
	return std::complex<RealType>(0.0, -1.0);
}

template <typename TargetParamsType,
          typename ModelType,
          typename WaveFunctionTransfType,
          typename LanczosSolverType,
          typename VectorWithOffsetType>
class TimeVectorsRungeKutta : public TimeVectorsBase<TargetParamsType,
                                                     ModelType,
                                                     WaveFunctionTransfType,
                                                     LanczosSolverType,
                                                     VectorWithOffsetType> {

	using BaseType                   = TimeVectorsBase<TargetParamsType,
	                                                   ModelType,
	                                                   WaveFunctionTransfType,
	                                                   LanczosSolverType,
	                                                   VectorWithOffsetType>;
	using PairType                   = typename BaseType::PairType;
	using RealType                   = typename TargetParamsType::RealType;
	using ModelHelperType            = typename ModelType::ModelHelperType;
	using LeftRightSuperType         = typename ModelHelperType::LeftRightSuperType;
	using BasisWithOperatorsType     = typename LeftRightSuperType::BasisWithOperatorsType;
	using SparseMatrixType           = typename BasisWithOperatorsType::SparseMatrixType;
	using ComplexOrRealType          = typename SparseMatrixType::value_type;
	using MatrixComplexOrRealType    = PsimagLite::Matrix<ComplexOrRealType>;
	using VectorComplexOrRealType    = typename PsimagLite::Vector<ComplexOrRealType>::Type;
	using VectorRealType             = typename PsimagLite::Vector<RealType>::Type;
	using TargetVectorType           = VectorComplexOrRealType;
	using VectorSizeType             = PsimagLite::Vector<SizeType>::Type;
	using VectorVectorWithOffsetType = typename PsimagLite::Vector<VectorWithOffsetType*>::Type;

public:

	TimeVectorsRungeKutta(const TargetParamsType&       tstStruct,
	                      VectorVectorWithOffsetType&   targetVectors,
	                      const ModelType&              model,
	                      const WaveFunctionTransfType& wft,
	                      const LeftRightSuperType&     lrs)
	    : BaseType(model, lrs, wft, "rungekutta")
	    , progress_("TimeVectorsRungeKutta")
	    , tstStruct_(tstStruct)
	    , targetVectors_(targetVectors)
	    , model_(model)
	    , wft_(wft)
	    , lrs_(lrs)
	{ }

	void calcTimeVectors(const VectorSizeType&       indices,
	                     RealType                    Eg,
	                     const VectorWithOffsetType& phi,
	                     const typename BaseType::ExtraData&) override
	{
		PsimagLite::OstringStream                     msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg << "using RungeKutta";

		RealType norma = norm(phi);
		if (norma < 1e-10)
			return;
		msg << " Norm of phi= " << norma;
		progress_.printline(msgg, std::cout);

		// set non-zero sectors
		for (SizeType i = 0; i < tstStruct_.times().size(); i++)
			*targetVectors_[i] = phi;

		for (SizeType ii = 0; ii < phi.sectors(); ii++) {
			SizeType i = phi.sector(ii);
			calcTimeVectors(indices, Eg, phi, i);
		}
	}

private:

	class FunctionForRungeKutta {

	public:

		FunctionForRungeKutta(const RealType&             E0,
		                      const RealType&             timeDirection,
		                      const LeftRightSuperType&   lrs,
		                      RealType                    currentTime,
		                      const ModelType&            model,
		                      const VectorWithOffsetType& phi,
		                      SizeType                    i0)
		    : E0_(E0)
		    , timeDirection_(timeDirection)
		    , p_(lrs.super().findPartitionNumber(phi.offset(i0)))
		    , aux_(p_, lrs)
		    , hc_(lrs,
		          ModelType::modelLinks(),
		          currentTime,
		          model.superOpHelper(),
		          model.ioIn())
		    , lanczosHelper_(model, hc_, aux_)
		{ }

		TargetVectorType operator()(const RealType&, const TargetVectorType& y) const
		{
			TargetVectorType x(y.size());
			lanczosHelper_.matrixVectorProduct(x, y);
			for (SizeType i = 0; i < x.size(); i++)
				x[i] -= E0_ * y[i];
			ComplexOrRealType icomplex
			    = minusOneOrMinusI(static_cast<ComplexOrRealType>(0));
			ComplexOrRealType tmp2 = timeDirection_ * icomplex;
			TargetVectorType  x2;
			x2 <= tmp2*       x;
			return x2;
		}

	private:

		RealType                                      E0_;
		RealType                                      timeDirection_;
		SizeType                                      p_;
		typename ModelHelperType::Aux                 aux_;
		typename ModelType::HamiltonianConnectionType hc_;
		typename LanczosSolverType::MatrixType        lanczosHelper_;
	}; // FunctionForRungeKutta

	void calcTimeVectors(const VectorSizeType&       indices,
	                     RealType                    Eg,
	                     const VectorWithOffsetType& phi,
	                     SizeType                    i0)
	{
		const VectorRealType& times = tstStruct_.times();
		SizeType              total = phi.effectiveSize(i0);
		TargetVectorType      phi0(total);
		phi.extract(phi0, i0);
		FunctionForRungeKutta f(
		    Eg, tstStruct_.timeDirection(), lrs_, this->time(), model_, phi, i0);

		RealType epsForRK = tstStruct_.tau() / (times.size() - 1.0);
		PsimagLite::RungeKutta<RealType, FunctionForRungeKutta, TargetVectorType>
		    rungeKutta(f, epsForRK);

		typename PsimagLite::Vector<TargetVectorType>::Type result;
		rungeKutta.solve(result, 0.0, times.size(), phi0);
		assert(result.size() == times.size());

		const SizeType n = indices.size();
		for (SizeType i = 0; i < n; ++i) {
			const SizeType ii = indices[i];
			assert(ii < targetVectors_.size());
			targetVectors_[ii]->setDataInSector(result[i], i0);
		}
	}

	PsimagLite::ProgressIndicator progress_;
	const TargetParamsType&       tstStruct_;
	VectorVectorWithOffsetType&   targetVectors_;
	const ModelType&              model_;
	const WaveFunctionTransfType& wft_;
	const LeftRightSuperType&     lrs_;
}; // class TimeVectorsRungeKutta
} // namespace Dmrg
/*@}*/
#endif
