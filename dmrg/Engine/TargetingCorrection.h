// SPDX-FileCopyrightText: Copyright (c) 2009-2013, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file TargetingCorrection.h
 *
 *  corrects the finite-size algorithm
 *  following PRB 72, 180403(R) (2005)
 *
 */

#ifndef CORRECTION_TARGETING_H
#define CORRECTION_TARGETING_H
#include "TargetParamsCorrection.h"
#include "TargetingBase.h"
#include <iostream>
#include <stdexcept>

namespace Dmrg {

template <typename LanczosSolverType_, typename VectorWithOffsetType_>
class TargetingCorrection : public TargetingBase<LanczosSolverType_, VectorWithOffsetType_> {

public:

	using LanczosSolverType      = LanczosSolverType_;
	using BaseType               = TargetingBase<LanczosSolverType, VectorWithOffsetType_>;
	using TargetingCommonType    = typename BaseType::TargetingCommonType;
	using MatrixVectorType       = typename BaseType::MatrixVectorType;
	using ModelType              = typename MatrixVectorType::ModelType;
	using RealType               = typename ModelType::RealType;
	using ModelHelperType        = typename ModelType::ModelHelperType;
	using LeftRightSuperType     = typename ModelHelperType::LeftRightSuperType;
	using BasisWithOperatorsType = typename LeftRightSuperType::BasisWithOperatorsType;
	using SparseMatrixType       = typename BasisWithOperatorsType::SparseMatrixType;
	using ComplexOrRealType      = typename SparseMatrixType::value_type;
	using TargetVectorType       = typename PsimagLite::Vector<ComplexOrRealType>::Type;
	using OperatorType           = typename BasisWithOperatorsType::OperatorType;
	using BasisType              = typename BasisWithOperatorsType::BasisType;
	using BlockType              = typename BasisType::BlockType;
	using WaveFunctionTransfType = typename BaseType::WaveFunctionTransfType;
	using VectorWithOffsetType   = typename WaveFunctionTransfType::VectorWithOffsetType;
	using TargetParamsType       = TargetParamsCorrection<ModelType>;
	using InputValidatorType     = typename ModelType::InputValidatorType;
	using VectorSizeType         = typename PsimagLite::Vector<SizeType>::Type;
	using QnType                 = typename BasisType::QnType;
	using StageEnumType          = typename TargetingCommonType::StageEnumType;
	using VectorRealType         = typename BaseType::VectorRealType;
	using CheckpointType         = typename BaseType::CheckpointType;

	TargetingCorrection(const LeftRightSuperType&     lrs,
	                    const CheckpointType&         checkPoint,
	                    const WaveFunctionTransfType& wft,
	                    const QnType&,
	                    InputValidatorType& io)
	    : BaseType(lrs, checkPoint, wft, 0)
	    , tstStruct_(io, "TargetingCorrection")
	    , gsWeight_(1 - tstStruct_.correctionA())
	    , progress_("TargetingCorrection")
	{ }

	SizeType sites() const override { return tstStruct_.sites(); }

	SizeType targets() const override { return 1; }

	RealType normSquared(SizeType i) const override
	{
		return PsimagLite::real(this->tv(i) * this->tv(i));
	}

	RealType weight(SizeType) const override
	{
		assert(this->common().aoe().noStageIs(StageEnumType::DISABLED));
		RealType gsWeight = 1;
		RealType weight1  = 0;
		computeAllWeights(gsWeight, weight1);
		return weight1;
	}

	RealType gsWeight() const override
	{
		RealType gsWeight = 1;
		RealType weight   = 0;
		computeAllWeights(gsWeight, weight);
		return gsWeight;
	}

	void evolve(const VectorRealType&,
	            ProgramGlobals::DirectionEnum direction,
	            const BlockType&              block1,
	            const BlockType&,
	            SizeType) override
	{
		if (direction == ProgramGlobals::DirectionEnum::INFINITE)
			return;

		this->common().setAllStagesTo(StageEnumType::WFT_NOADVANCE);
		this->common().computeCorrection(direction, block1);
		bool doBorderIfBorder = false;
		this->common().cocoon(block1, direction, doBorderIfBorder);
	}

	void read(typename TargetingCommonType::IoInputType& io, PsimagLite::String prefix) override
	{
		this->common().read(io, prefix);
	}

	void write(const typename PsimagLite::Vector<SizeType>::Type& block,
	           PsimagLite::IoSelector::Out&                       io,
	           PsimagLite::String                                 prefix) const override
	{
		this->common().write(io, block, prefix);
	}

private:

	void computeAllWeights(RealType& gsWeight, RealType& weight) const
	{
		assert(1 == this->common().aoe().tvs());

		RealType norma = norm(this->tv(0));
		weight         = (norma > 1e-6) ? tstStruct_.correctionA() : 0;

		gsWeight = 1 - weight;

		if (gsWeight >= gsWeight_)
			return; // <--- EARLY EXIT HERE

		assert(weight > 1e-6);
		RealType factor = (1 - gsWeight_) / weight;
		weight *= factor;

		gsWeight = gsWeight_;
	}

	TargetParamsType              tstStruct_;
	const RealType                gsWeight_;
	PsimagLite::ProgressIndicator progress_;
}; // class TargetingCorrection
} // namespace Dmrg
/*@}*/
#endif
