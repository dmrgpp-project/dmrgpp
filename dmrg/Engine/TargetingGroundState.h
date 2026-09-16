// SPDX-FileCopyrightText: Copyright (c) 2009, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 2.0.0]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file TargetingGroundState.h
 *
 *  targets the ground state
 *
 */

#ifndef TARGETING_GS_H
#define TARGETING_GS_H
#include "ApplyOperatorLocal.h"
#include "TargetParamsGroundState.h"
#include "TargetingBase.h"
#include <PsimagLite/ParametersForSolver.h>
#include <PsimagLite/PsimagLite.h>
#include <iostream>
#include <stdexcept>

namespace Dmrg {

template <typename LanczosSolverType_, typename VectorWithOffsetType_>
class TargetingGroundState : public TargetingBase<LanczosSolverType_, VectorWithOffsetType_> {

public:

	using LanczosSolverType      = LanczosSolverType_;
	using BaseType               = TargetingBase<LanczosSolverType, VectorWithOffsetType_>;
	using TargetingCommonType    = typename BaseType::TargetingCommonType;
	using MatrixVectorType       = typename BaseType::MatrixVectorType;
	using CheckpointType         = typename BaseType::CheckpointType;
	using ModelType              = typename MatrixVectorType::ModelType;
	using RealType               = typename ModelType::RealType;
	using ModelHelperType        = typename ModelType::ModelHelperType;
	using LeftRightSuperType     = typename ModelHelperType::LeftRightSuperType;
	using BasisWithOperatorsType = typename LeftRightSuperType::BasisWithOperatorsType;
	using SparseMatrixType       = typename BasisWithOperatorsType::SparseMatrixType;
	using OperatorType           = typename BasisWithOperatorsType::OperatorType;
	using BasisType              = typename BasisWithOperatorsType::BasisType;
	using ComplexOrRealType      = typename SparseMatrixType::value_type;
	using VectorType             = typename PsimagLite::Vector<ComplexOrRealType>::Type;
	using BlockType              = typename BasisType::BlockType;
	using WaveFunctionTransfType = typename BaseType::WaveFunctionTransfType;
	using VectorWithOffsetType   = typename WaveFunctionTransfType::VectorWithOffsetType;
	using TargetVectorType       = VectorType;
	using TargetParamsType       = TargetParamsGroundState<ModelType>;
	using InputValidatorType     = typename ModelType::InputValidatorType;
	using VectorSizeType         = typename PsimagLite::Vector<SizeType>::Type;
	using QnType                 = typename BasisType::QnType;
	using VectorRealType         = typename BaseType::VectorRealType;
	using VectorVectorVectorWithOffsetType =
	    typename TargetingCommonType::VectorVectorVectorWithOffsetType;
	using VectorVectorWithOffsetType = typename PsimagLite::Vector<VectorWithOffsetType*>::Type;

	TargetingGroundState(const LeftRightSuperType&     lrs,
	                     const CheckpointType&         checkPoint,
	                     const WaveFunctionTransfType& wft,
	                     const QnType&,
	                     InputValidatorType&)
	    : BaseType(lrs, checkPoint, wft, 0)
	    , tstStruct_("TargetingGroundState")
	    , progress_("TargetingGroundState")
	{ }

	SizeType sites() const override { return tstStruct_.sites(); }

	SizeType targets() const override { return 0; }

	RealType weight(SizeType) const override
	{
		throw PsimagLite::RuntimeError("GST: What are you doing here?\n");
	}

	RealType gsWeight() const override
	{
		const VectorVectorVectorWithOffsetType& v = this->common().aoe().psiConst();
		const SizeType                          n = v.size();

		RealType sum = 0;
		for (SizeType i = 0; i < n; ++i) {
			const VectorVectorWithOffsetType& vv = v[i];
			const SizeType                    m  = vv.size();
			for (SizeType j = 0; j < m; ++j) {
				const VectorWithOffsetType* vvv = vv[j];
				if (!vvv)
					continue;
				if (norm(*vvv) > 1e-3)
					++sum;
			}
		}

		return 1.0 / sum;
	}

	SizeType size() const override { return 0; }

	void evolve(const VectorRealType&,
	            ProgramGlobals::DirectionEnum direction,
	            const BlockType&              block1,
	            const BlockType&,
	            SizeType) override
	{
		bool doBorderIfBorder = true;
		this->common().cocoon(block1, direction, doBorderIfBorder);
	}

	void write(const typename PsimagLite::Vector<SizeType>::Type& block,
	           PsimagLite::IoSelector::Out&                       io,
	           PsimagLite::String                                 prefix) const override
	{
		this->common().write(io, block, prefix);
	}

	void read(typename TargetingCommonType::IoInputType& io, PsimagLite::String prefix) override
	{
		this->common().readGSandNGSTs(io, prefix, "GroundState");
	}

private:

	TargetParamsType              tstStruct_;
	PsimagLite::ProgressIndicator progress_;

}; // class TargetingGroundState
} // namespace Dmrg
/*@}*/
#endif
