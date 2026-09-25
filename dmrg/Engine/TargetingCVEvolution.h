// SPDX-FileCopyrightText: Copyright (c) 2009-2014-2021, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 6.]
// [DMRG++ authors; see AUTHORS.]

#ifndef TARGETING_CVEvolution_H
#define TARGETING_CVEvolution_H

#include "BlockDiagonalMatrix.h"
#include "CorrectionVectorSkeleton.h"
#include "ProgramGlobals.h"
#include "TargetParamsCorrectionVector.h"
#include "TargetingBase.h"
#include "TimeVectorsKrylov.h"
#include <PsimagLite/ParametersForSolver.h>
#include <PsimagLite/PredicateAwesome.h>
#include <PsimagLite/ProgressIndicator.h>
#include <iostream>

namespace Dmrg {

template <typename LanczosSolverType_, typename VectorWithOffsetType_>
class TargetingCVEvolution : public TargetingBase<LanczosSolverType_, VectorWithOffsetType_> {

	enum
	{
		BORDER_NEITHER,
		BORDER_LEFT,
		BORDER_RIGHT
	};

public:

	using LanczosSolverType      = LanczosSolverType_;
	using BaseType               = TargetingBase<LanczosSolverType, VectorWithOffsetType_>;
	using TargetingCommonType    = typename BaseType::TargetingCommonType;
	using PairType               = std::pair<SizeType, SizeType>;
	using OptionsType            = typename BaseType::OptionsType;
	using MatrixVectorType       = typename BaseType::MatrixVectorType;
	using CheckpointType         = typename BaseType::CheckpointType;
	using ModelType              = typename MatrixVectorType::ModelType;
	using RealType               = typename ModelType::RealType;
	using OperatorsType          = typename ModelType::OperatorsType;
	using ModelHelperType        = typename ModelType::ModelHelperType;
	using LeftRightSuperType     = typename ModelHelperType::LeftRightSuperType;
	using BasisWithOperatorsType = typename LeftRightSuperType::BasisWithOperatorsType;
	using VectorSizeType         = PsimagLite::Vector<SizeType>::Type;
	using WaveFunctionTransfType = typename BaseType::WaveFunctionTransfType;
	using VectorWithOffsetType   = typename WaveFunctionTransfType::VectorWithOffsetType;
	using ComplexOrRealType      = typename VectorWithOffsetType::value_type;
	using TargetVectorType       = typename VectorWithOffsetType::VectorType;
	using VectorRealType         = typename PsimagLite::Vector<RealType>::Type;
	using OperatorType           = typename BasisWithOperatorsType::OperatorType;
	using BasisType              = typename BasisWithOperatorsType::BasisType;
	using TargetParamsType       = TargetParamsCorrectionVector<ModelType>;
	using BlockType              = typename BasisType::BlockType;
	using TimeSerializerType     = typename TargetingCommonType::TimeSerializerType;
	using SparseMatrixType       = typename OperatorType::StorageType;
	using InputValidatorType     = typename ModelType::InputValidatorType;
	using QnType                 = typename BasisType::QnType;
	using StageEnumType          = typename TargetingCommonType::StageEnumType;
	using CorrectionVectorSkeletonType = CorrectionVectorSkeleton<LanczosSolverType,
	                                                              VectorWithOffsetType,
	                                                              BaseType,
	                                                              TargetParamsType>;
	typedef
	    typename TargetingCommonType::ApplyOperatorExpressionType ApplyOperatorExpressionType;
	using TimeVectorsBaseType = typename ApplyOperatorExpressionType::TimeVectorsBaseType;

	TargetingCVEvolution(const LeftRightSuperType&     lrs,
	                     const CheckpointType&         checkPoint,
	                     const WaveFunctionTransfType& wft,
	                     const QnType&,
	                     InputValidatorType& ioIn)
	    : BaseType(lrs, checkPoint, wft, 0)
	    , tstStruct_(ioIn, "TargetingCVEvolution", checkPoint.model())
	    , wft_(wft)
	    , progress_("TargetingCVEvolution")
	    , counter_(0)
	    , almostDone_(0)
	    , skeleton_(ioIn, tstStruct_, checkPoint.model(), lrs, this->common().aoe().energy())
	    , weight_(targets())
	    , gsWeight_(tstStruct_.gsWeight())
	{
		if (!wft.isEnabled())
			err("TST needs an enabled wft\n");

		if (tstStruct_.sites() == 0)
			err("TST needs at least one TSPSite\n");

		if (gsWeight_ < 0 || gsWeight_ >= 1)
			err("gsWeight_ must be in [0, 1)\n");

		SizeType n      = weight_.size();
		RealType factor = (1.0 - gsWeight_) / n;
		for (SizeType i = 0; i < n; ++i)
			weight_[i] = factor;
	}

	SizeType sites() const override { return tstStruct_.sites(); }

	SizeType targets() const override
	{
		return PsimagLite::IsComplexNumber<ComplexOrRealType>::True ? 3 : 5;
	}

	RealType weight(SizeType i) const override
	{
		assert(!this->common().aoe().allStages(StageEnumType::DISABLED));
		return weight_[i];
	}

	RealType gsWeight() const override
	{
		if (this->common().aoe().allStages(StageEnumType::DISABLED))
			return 1.0;
		return gsWeight_;
	}

	bool includeGroundStage() const override
	{
		if (!this->common().aoe().noStageIs(StageEnumType::DISABLED))
			return true;
		bool b = (fabs(gsWeight_) > 1e-6);
		return b;
	}

	void evolve(const VectorRealType&         energies,
	            ProgramGlobals::DirectionEnum direction,
	            const BlockType&              block1,
	            const BlockType&,
	            SizeType loopNumber) override
	{
		assert(block1.size() > 0);
		SizeType site = block1[0];
		assert(energies.size() > 0);
		RealType Eg = energies[0];
		evolveInternal(Eg, direction, block1, loopNumber);
		SizeType numberOfSites = this->lrs().super().block().size();

		if (site > 1 && site < numberOfSites - 2)
			return;

		if (direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) {
			if (site == 1)
				return;
		} else {
			if (site == numberOfSites - 2)
				return;
		}

		SizeType  x = (site == 1) ? 0 : numberOfSites - 1;
		BlockType block(1, x);
		evolveInternal(Eg, direction, block, loopNumber);
	}

	bool end() const override { return (almostDone_ >= 2); }

	void read(typename TargetingCommonType::IoInputType& io, PsimagLite::String prefix) override
	{
		this->common().readGSandNGSTs(io, prefix, "CVEvolution");
	}

	void write(const VectorSizeType&        block,
	           PsimagLite::IoSelector::Out& io,
	           PsimagLite::String           prefix) const override
	{
		PsimagLite::OstringStream                     msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg << "Saving state...";
		progress_.printline(msgg, std::cout);

		this->common().write(io, block, prefix);
		this->common().writeNGSTs(io, prefix, block, "TimeStep");
	}

private:

	void evolveInternal(RealType                      Eg,
	                    ProgramGlobals::DirectionEnum direction,
	                    const BlockType&              block1,
	                    SizeType                      loopNumber)
	{
		if (direction == ProgramGlobals::DirectionEnum::INFINITE)
			return;
		VectorWithOffsetType phiNew;
		assert(block1.size() > 0);
		SizeType site = block1[0];
		this->common().aoeNonConst().getPhi(
		    &phiNew, Eg, direction, site, loopNumber, tstStruct_);

		if (phiNew.size() == 0)
			return;

		this->tvNonConst(0) = phiNew;
		VectorWithOffsetType bogusTv;

		const SizeType currentTimeStep
		    = this->common().aoe().timeVectors().currentTimeStep();
		if (currentTimeStep == 0) {
			if (PsimagLite::IsComplexNumber<ComplexOrRealType>::True) {
				skeleton_.calcDynVectors(phiNew, this->tvNonConst(1), bogusTv);
				skeleton_.calcDynVectors(this->tv(1), this->tvNonConst(2), bogusTv);
			} else {

				skeleton_.calcDynVectors(
				    phiNew, this->tvNonConst(1), this->tvNonConst(2));

				skeleton_.calcDynVectors(this->tv(1),
				                         this->tv(2),
				                         this->tvNonConst(3),
				                         this->tvNonConst(4));
			}
		} else {
			bool timeHasAdvanced = (counter_ != currentTimeStep
			                        && currentTimeStep < tstStruct_.nForFraction());

			if (counter_ != currentTimeStep && !timeHasAdvanced) {
				std::cout << __FILE__ << " is now DONE\n";
				std::cerr << __FILE__ << " is now DONE\n";
				++almostDone_;
			}

			if (PsimagLite::IsComplexNumber<ComplexOrRealType>::True) {

				const SizeType advanceIndex = (timeHasAdvanced) ? 2 : 1;
				// wft tv1
				this->common().aoe().wftOneVector(
				    bogusTv, this->tv(advanceIndex), site);
				this->tvNonConst(1) = bogusTv;
				skeleton_.calcDynVectors(this->tv(1), this->tvNonConst(2), bogusTv);
			} else {

				VectorWithOffsetType bogusTv2;
				const SizeType       advanceIndex = (timeHasAdvanced) ? 3 : 1;

				this->common().aoe().wftOneVector(
				    bogusTv, this->tv(advanceIndex), site);
				const SizeType advanceIndexp1 = advanceIndex + 1;
				this->common().aoe().wftOneVector(
				    bogusTv2, this->tv(advanceIndexp1), site);
				this->tvNonConst(1) = bogusTv;
				this->tvNonConst(2) = bogusTv2;

				skeleton_.calcDynVectors(this->tv(1),
				                         this->tv(2),
				                         this->tvNonConst(3),
				                         this->tvNonConst(4));
			}
		}

		counter_ = currentTimeStep;
		auto* ptr
		    = const_cast<TimeVectorsBaseType*>(&this->common().aoeNonConst().timeVectors());
		ptr->setCurrentTime(counter_);

		bool doBorderIfBorder = true;
		this->common().cocoon(block1, direction, doBorderIfBorder);

		this->common().printNormsAndWeights(gsWeight_, weight_);
	}

	TargetParamsType              tstStruct_;
	const WaveFunctionTransfType& wft_;
	PsimagLite::ProgressIndicator progress_;
	SizeType                      counter_;
	SizeType                      almostDone_;
	CorrectionVectorSkeletonType  skeleton_;
	VectorRealType                weight_;
	RealType                      gsWeight_;
}; // class TargetingCVEvolution
} // namespace Dmrg

#endif
