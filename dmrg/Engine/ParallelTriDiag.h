// SPDX-FileCopyrightText: Copyright (c) 2009,-2014 UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/
/** \file ParallelTriDiag.h
 */

#ifndef PARALLEL_TRIDIAG_H
#define PARALLEL_TRIDIAG_H

#include <PsimagLite/Concurrency.h>
#include <PsimagLite/NotMpi.h>

namespace Dmrg {

template <typename ModelType, typename LanczosSolverType, typename VectorWithOffsetType>
class ParallelTriDiag {

	using ModelHelperType        = typename ModelType::ModelHelperType;
	using LeftRightSuperType     = typename ModelHelperType::LeftRightSuperType;
	using BasisWithOperatorsType = typename LeftRightSuperType::BasisWithOperatorsType;
	using SparseMatrixType       = typename BasisWithOperatorsType::SparseMatrixType;
	using ComplexOrRealType      = typename SparseMatrixType::value_type;
	using RealType               = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using TridiagonalMatrixType  = typename LanczosSolverType::TridiagonalMatrixType;
	using InputValidatorType     = typename ModelType::InputValidatorType;
	using ConcurrencyType        = PsimagLite::Concurrency;

public:

	using TargetVectorType        = typename PsimagLite::Vector<ComplexOrRealType>::Type;
	using MatrixComplexOrRealType = PsimagLite::Matrix<ComplexOrRealType>;
	using VectorMatrixFieldType   = typename PsimagLite::Vector<MatrixComplexOrRealType>::Type;

	ParallelTriDiag(const VectorWithOffsetType&                  phi,
	                VectorMatrixFieldType&                       T,
	                VectorMatrixFieldType&                       V,
	                typename PsimagLite::Vector<SizeType>::Type& steps,
	                const LeftRightSuperType&                    lrs,
	                RealType                                     currentTime,
	                const ModelType&                             model,
	                InputValidatorType&                          io)
	    : phi_(phi)
	    , T_(T)
	    , V_(V)
	    , steps_(steps)
	    , lrs_(lrs)
	    , currentTime_(currentTime)
	    , model_(model)
	    , io_(io)
	{ }

	SizeType tasks() const { return phi_.sectors(); }

	void doTask(SizeType ii, SizeType)
	{
		SizeType i = phi_.sector(ii);
		steps_[ii] = triDiag(phi_, T_[ii], V_[ii], i);
	}

private:

	SizeType triDiag(const VectorWithOffsetType& phi,
	                 MatrixComplexOrRealType&    T,
	                 MatrixComplexOrRealType&    V,
	                 SizeType                    i0)
	{
		const SizeType                p = lrs_.super().findPartitionNumber(phi.offset(i0));
		typename ModelHelperType::Aux aux(p, lrs_);
		typename ModelType::HamiltonianConnectionType hc(lrs_,
		                                                 ModelType::modelLinks(),
		                                                 currentTime_,
		                                                 model_.superOpHelper(),
		                                                 model_.ioIn());
		typename LanczosSolverType::MatrixType        lanczosHelper(model_, hc, aux);

		typename LanczosSolverType::ParametersSolverType params(io_, "Tridiag");
		params.lotaMemory = true;

		LanczosSolverType lanczosSolver(lanczosHelper, params);

		TridiagonalMatrixType ab;
		SizeType              total = phi.effectiveSize(i0);
		TargetVectorType      phi2(total);
		phi.extract(phi2, i0);
		lanczosSolver.decomposition(phi2, ab);
		ab.buildDenseMatrix(T);

		lanczosSolver.lanczosVectorsSwap(V);

		return lanczosSolver.steps();
	}

	const VectorWithOffsetType&                  phi_;
	VectorMatrixFieldType&                       T_;
	VectorMatrixFieldType&                       V_;
	typename PsimagLite::Vector<SizeType>::Type& steps_;
	const LeftRightSuperType&                    lrs_;
	RealType                                     currentTime_;
	const ModelType&                             model_;
	InputValidatorType&                          io_;
}; // class ParallelTriDiag
} // namespace Dmrg

/*@}*/
#endif // PARALLEL_TRIDIAG_H
