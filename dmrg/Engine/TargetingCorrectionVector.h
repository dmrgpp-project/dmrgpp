// SPDX-FileCopyrightText: Copyright (c) 2009-2011, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 2.0.0]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file TargetingCorrectionVector.h
 *
 * Implements the targeting required by
 * the correction targeting method
 *
 */

#ifndef TARGETING_CORRECTION_VECTOR_H
#define TARGETING_CORRECTION_VECTOR_H

#include "CorrectionVectorFunction.h"
#include "CorrectionVectorSkeleton.h"
#include "ParallelTriDiag.h"
#include "TargetParamsCorrectionVector.h"
#include "TargetingBase.h"
#include "VectorWithOffsets.h"
#include <PsimagLite/FreqEnum.h>
#include <PsimagLite/NoPthreadsNg.h>
#include <PsimagLite/ParametersForSolver.h>
#include <PsimagLite/ProgressIndicator.h>

namespace Dmrg {

template <typename LanczosSolverType_, typename VectorWithOffsetType_>
class TargetingCorrectionVector : public TargetingBase<LanczosSolverType_, VectorWithOffsetType_> {

	using LanczosSolverType = LanczosSolverType_;
	using BaseType          = TargetingBase<LanczosSolverType, VectorWithOffsetType_>;

public:

	using TargetingCommonType    = typename BaseType::TargetingCommonType;
	using MatrixVectorType       = typename BaseType::MatrixVectorType;
	using ModelType              = typename MatrixVectorType::ModelType;
	using RealType               = typename ModelType::RealType;
	using OptionsType            = typename BaseType::OptionsType;
	using CheckpointType         = typename BaseType::CheckpointType;
	using VectorRealType         = typename PsimagLite::Vector<RealType>::Type;
	using OperatorsType          = typename ModelType::OperatorsType;
	using ModelHelperType        = typename ModelType::ModelHelperType;
	using LeftRightSuperType     = typename ModelHelperType::LeftRightSuperType;
	using BasisWithOperatorsType = typename LeftRightSuperType::BasisWithOperatorsType;
	using OperatorType           = typename BasisWithOperatorsType::OperatorType;
	using BasisType              = typename BasisWithOperatorsType::BasisType;
	using SparseMatrixType       = typename BasisWithOperatorsType::SparseMatrixType;
	using ComplexOrRealType      = typename SparseMatrixType::value_type;
	using TargetParamsType       = TargetParamsCorrectionVector<ModelType>;
	using BlockType              = typename BasisType::BlockType;
	using WaveFunctionTransfType = typename BaseType::WaveFunctionTransfType;
	using VectorWithOffsetType   = typename WaveFunctionTransfType::VectorWithOffsetType;
	using VectorType             = typename VectorWithOffsetType::VectorType;
	using TargetVectorType       = VectorType;
	using TimeSerializerType     = typename TargetingCommonType::TimeSerializerType;
	using TridiagonalMatrixType  = typename LanczosSolverType::TridiagonalMatrixType;
	using DenseMatrixType        = PsimagLite::Matrix<typename VectorType::value_type>;
	using DenseMatrixRealType    = PsimagLite::Matrix<RealType>;
	using PostProcType           = typename LanczosSolverType::PostProcType;
	using LanczosMatrixType      = typename LanczosSolverType::MatrixType;
	using CorrectionVectorFunctionType
	    = CorrectionVectorFunction<LanczosMatrixType, TargetParamsType>;
	using ParallelTriDiagType
	    = ParallelTriDiag<ModelType, LanczosSolverType, VectorWithOffsetType>;
	using MatrixComplexOrRealType      = typename ParallelTriDiagType::MatrixComplexOrRealType;
	using VectorMatrixFieldType        = typename ParallelTriDiagType::VectorMatrixFieldType;
	using VectorSizeType               = typename PsimagLite::Vector<SizeType>::Type;
	using VectorVectorRealType         = typename PsimagLite::Vector<VectorRealType>::Type;
	using InputValidatorType           = typename ModelType::InputValidatorType;
	using CorrectionVectorSkeletonType = CorrectionVectorSkeleton<LanczosSolverType,
	                                                              VectorWithOffsetType,
	                                                              BaseType,
	                                                              TargetParamsType>;
	using QnType                       = typename BasisType::QnType;

	TargetingCorrectionVector(const LeftRightSuperType&     lrs,
	                          const CheckpointType&         checkPoint,
	                          const WaveFunctionTransfType& wft,
	                          const QnType&,
	                          InputValidatorType& ioIn)
	    : BaseType(lrs, checkPoint, wft, 1)
	    , tstStruct_(ioIn, "TargetingCorrectionVector", checkPoint.model())
	    , ioIn_(ioIn)
	    , progress_("TargetingCorrectionVector")
	    , gsWeight_(1.0)
	    , correctionEnabled_(false)
	    , skeleton_(ioIn_, tstStruct_, checkPoint.model(), lrs, this->common().aoe().energy())
	{
		if (!wft.isEnabled())
			err("TargetingCorrectionVector needs wft\n");
	}

	SizeType sites() const override { return tstStruct_.sites(); }

	SizeType targets() const override { return 4; }

	RealType weight(SizeType i) const override
	{
		assert(i < weight_.size());
		return weight_[i];
	}

	RealType gsWeight() const override
	{
		if (!correctionEnabled_)
			return 1.0;
		return gsWeight_;
	}

	SizeType size() const override
	{
		if (!correctionEnabled_)
			return 0;
		return BaseType::size();
	}

	void evolve(const VectorRealType&         energies,
	            ProgramGlobals::DirectionEnum direction,
	            const BlockType&              block1,
	            const BlockType&              block2,
	            SizeType                      loopNumber) override
	{
		if (block1.size() != 1 || block2.size() != 1) {
			PsimagLite::String str(__FILE__);
			str += " " + ttos(__LINE__) + "\n";
			str += "evolve only blocks of one site supported\n";
			throw PsimagLite::RuntimeError(str.c_str());
		}

		assert(energies.size() > 0);
		RealType Eg   = energies[0];
		SizeType site = block1[0];
		evolve(Eg, direction, site, loopNumber);

		this->common().printNormsAndWeights(gsWeight_, weight_);

		// corner case
		SizeType numberOfSites = this->lrs().super().block().size();
		SizeType site2         = numberOfSites;

		if (site == 1 && direction == ProgramGlobals::DirectionEnum::EXPAND_ENVIRON)
			site2 = 0;
		if (site == numberOfSites - 2
		    && direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM)
			site2 = numberOfSites - 1;
		if (site2 == numberOfSites)
			return;
		evolve(Eg, direction, site2, loopNumber);
	}

	void write(const typename PsimagLite::Vector<SizeType>::Type& block,
	           PsimagLite::IoSelector::Out&                       io,
	           PsimagLite::String                                 prefix) const override
	{
		this->common().write(io, block, prefix);
		this->common().writeNGSTs(io, prefix, block, "CorrectionVector");
	}

	void read(typename TargetingCommonType::IoInputType& io, PsimagLite::String prefix) override
	{
		this->common().readGSandNGSTs(io, prefix, "CorrectionVector");
		setWeights();
	}

private:

	void evolve(RealType                      Eg,
	            ProgramGlobals::DirectionEnum direction,
	            SizeType                      site,
	            SizeType                      loopNumber)
	{
		VectorWithOffsetType phiNew;
		SizeType             count = this->common().aoeNonConst().getPhi(
                    &phiNew, Eg, direction, site, loopNumber, tstStruct_);

		if (direction != ProgramGlobals::DirectionEnum::INFINITE) {
			correctionEnabled_ = true;
			typename PsimagLite::Vector<SizeType>::Type block1(1, site);
			addCorrection(direction, block1);
		}

		if (count == 0)
			return;

		this->tvNonConst(1) = phiNew;
		skeleton_.calcDynVectors(this->tv(1), this->tvNonConst(2), this->tvNonConst(3));

		setWeights();

		VectorSizeType block(1, site);
		bool           doBorderIfBorder = false;
		this->common().cocoon(block, direction, doBorderIfBorder);
	}

	void setWeights()
	{
		gsWeight_ = tstStruct_.gsWeight();

		RealType sum = 0;
		weight_.resize(this->common().aoe().tvs());

		for (SizeType r = 1; r < weight_.size(); r++) {
			weight_[r] = 1;
			sum += weight_[r];
		}
		for (SizeType r = 0; r < weight_.size(); r++)
			weight_[r] *= (1.0 - gsWeight_) / sum;

		if (PsimagLite::IsComplexNumber<ComplexOrRealType>::True) {
			sum = 0;
			for (SizeType r = 1; r < weight_.size(); r++) {
				weight_[r] = 1;
				sum += weight_[r];
			}
			weight_[weight_.size() - 1] = 0;
			for (SizeType r = 1; r < weight_.size(); r++)
				weight_[r] *= (1.0 - gsWeight_) * r / sum;
		}
	}

	RealType dynWeightOf(VectorType& v, const VectorType& w) const
	{
		RealType sum = 0;
		for (SizeType i = 0; i < v.size(); i++) {
			RealType tmp = PsimagLite::real(v[i] * w[i]);
			sum += tmp * tmp;
		}
		return sum;
	}

	void addCorrection(ProgramGlobals::DirectionEnum direction, const BlockType& block1)
	{
		if (tstStruct_.correctionA() == 0)
			return;
		weight_.resize(1);
		weight_[0] = tstStruct_.correctionA();
		this->common().computeCorrection(direction, block1);
		gsWeight_ = 1.0 - weight_[0];
	}

	TargetParamsType                            tstStruct_;
	InputValidatorType&                         ioIn_;
	PsimagLite::ProgressIndicator               progress_;
	RealType                                    gsWeight_;
	bool                                        correctionEnabled_;
	typename PsimagLite::Vector<RealType>::Type weight_;
	CorrectionVectorSkeletonType                skeleton_;
}; // class TargetingCorrectionVector
} // namespace
/*@}*/
#endif // TARGETING_CORRECTION_VECTOR_H
