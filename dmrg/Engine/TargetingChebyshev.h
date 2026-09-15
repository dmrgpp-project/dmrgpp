// SPDX-FileCopyrightText: Copyright (c) 2009-2014, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [by G.A., Oak Ridge National Laboratory]


#ifndef TARGETING_CHEBYSHEV_H
#define TARGETING_CHEBYSHEV_H

#include "BlockDiagonalMatrix.h"
#include "OracleChebyshev.h"
#include "ProgramGlobals.h"
#include "TargetParamsTimeStep.h"
#include "TargetingBase.h"
#include "TimeVectorsChebyshev.h"
#include <PsimagLite/ParametersForSolver.h>
#include <PsimagLite/ProgressIndicator.h>
#include <iostream>

namespace Dmrg {

template <typename LanczosSolverType_, typename VectorWithOffsetType_>
class TargetingChebyshev : public TargetingBase<LanczosSolverType_, VectorWithOffsetType_> {

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
	using OptionsType            = typename BaseType::OptionsType;
	using PairType               = std::pair<SizeType, SizeType>;
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
	using TargetParamsType       = TargetParamsTimeStep<ModelType>;
	using BlockType              = typename BasisType::BlockType;
	using TimeSerializerType     = typename TargetingCommonType::TimeSerializerType;
	using SparseMatrixType       = typename OperatorType::StorageType;
	using InputValidatorType     = typename ModelType::InputValidatorType;
	using QnType                 = typename BasisType::QnType;
	typedef
	    typename TargetingCommonType::ApplyOperatorExpressionType ApplyOperatorExpressionType;
	using ApplyOperatorType = typename ApplyOperatorExpressionType::ApplyOperatorType;
	using StageEnumType     = typename TargetingCommonType::StageEnumType;

	TargetingChebyshev(const LeftRightSuperType&     lrs,
	                   const CheckpointType&         checkPoint,
	                   const WaveFunctionTransfType& wft,
	                   const QnType&,
	                   InputValidatorType& ioIn)
	    : BaseType(lrs, checkPoint, wft, 0)
	    , tstStruct_(ioIn, "TargetingChebyshev", checkPoint.model())
	    , wft_(wft)
	    , progress_("TargetingChebyshev")
	    , weight_(tstStruct_.times().size())
	    , tvEnergy_(weight_.size(), 0.0)
	    , gsWeight_(tstStruct_.gsWeight())
	{
		if (!wft.isEnabled())
			err("TST needs an enabled wft\n");
		if (tstStruct_.sites() == 0)
			err("TST needs at least one TSPSite\n");

		SizeType nops    = tstStruct_.sites();
		SizeType linSize = checkPoint.model().superGeometry().numberOfSites();
		for (SizeType i = 0; i < nops; ++i)
			if (tstStruct_.sites(i) == 0 || tstStruct_.sites(i) == linSize - 1)
				err("TargetingChebyshev: FATAL: No application of operators at "
				    "borders\n");

		RealType tau = tstStruct_.tau();
		RealType sum = 0;
		SizeType n   = tstStruct_.times().size();

		if (n < 3)
			throw PsimagLite::RuntimeError(
			    "At least 3 Chebyshev vectors need to be targets\n");

		RealType factor = (n + 4.0) / (n + 2.0);
		factor *= (1.0 - gsWeight_);
		for (SizeType i = 0; i < n; i++) {
			tstStruct_.times()[i] = i * tau / (n - 1);
			weight_[i]            = factor / (n + 4);
			sum += weight_[i];
		}

		sum -= weight_[0];
		sum -= weight_[n - 1];
		weight_[0] = weight_[n - 1] = 2 * factor / (n + 4);
		sum += weight_[n - 1];
		sum += weight_[0];

		gsWeight_ = 1.0 - sum;
		sum += gsWeight_;
		assert(fabs(sum - 1.0) < 1e-5);

		this->common().aoeNonConst().initTimeVectors(tstStruct_, ioIn);
	}

	SizeType sites() const override { return tstStruct_.sites(); }

	SizeType targets() const override { return tstStruct_.times().size(); }

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
		assert(energies.size() > 0);
		RealType Eg = energies[0];
		evolveInternal(Eg, direction, block1, loopNumber);
		bool doBorderIfBorder = true;
		this->common().cocoon(block1, direction, doBorderIfBorder); // in-situ
	}

	bool end() const override
	{
		return (tstStruct_.maxTime() != 0
		        && this->common().aoe().timeVectors().time() >= tstStruct_.maxTime());
	}

	void read(typename TargetingCommonType::IoInputType& io, PsimagLite::String prefix) override
	{
		this->common().readGSandNGSTs(io, prefix, "Chebyshev");
	}

	void write(const VectorSizeType&        block,
	           PsimagLite::IoSelector::Out& io,
	           PsimagLite::String           prefix) const override
	{
		PsimagLite::OstringStream                     msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg << "Saving state...";
		progress_.printline(msgg, std::cout);

		assert(block.size() > 0);

		this->common().write(io, block, prefix);
		this->common().writeNGSTs(io, prefix, block, "Chebyshev");
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
		this->common().aoeNonConst().getPhi(
		    &phiNew, Eg, direction, block1[0], loopNumber, tstStruct_);

		if (phiNew.size() == 0)
			return;

		VectorSizeType indices(tstStruct_.times().size());
		for (SizeType i = 0; i < tstStruct_.times().size(); ++i)
			indices[i] = i;

		bool allOperatorsApplied
		    = (this->common().aoe().noStageIs(StageEnumType::DISABLED)
		       && this->common().aoe().noStageIs(StageEnumType::OPERATOR));

		assert(0 < block1.size());

		bool isLastCall = true;
		this->common().aoeNonConst().calcTimeVectors(
		    indices,
		    Eg,
		    phiNew,
		    direction,
		    allOperatorsApplied,
		    false, // don't wft or advance indices[0]
		    block1,
		    isLastCall);

		assert(phiNew.offset(0) == this->tv(1).offset(0));

		const OptionsType& options              = this->model().params().options;
		const bool         normalizeTimeVectors = (options.isSet("normalizeVectors")
                                                   && !options.isSet("neverNormalizeVectors"));

		assert(phiNew.offset(0) == this->tv(1).offset(0));

		if (normalizeTimeVectors)
			this->common().normalizeTimeVectors();

		this->common().printNormsAndWeights(gsWeight_, weight_);

		assert(phiNew.offset(0) == this->tv(1).offset(0));
	}

	void oracleChebyshev(SizeType site, SizeType systemOrEviron) const
	{
		using OracleChebyshevType = OracleChebyshev<TargetingCommonType, TargetParamsType>;
		OracleChebyshevType oracle(BaseType::model(),
		                           BaseType::lrs(),
		                           this->common().aoe().currentTime(),
		                           tstStruct_,
		                           this->common().aoe().energy());

		OperatorType A = BaseType::model().naturalOperator("c", 0, 0);
		oracle(3, this->common(), systemOrEviron, site, A, ApplyOperatorType::BORDER_NO);
	}

	void printChebyshev() const
	{
		for (SizeType i = 0; i < this->tv().size(); i++)
			printChebyshev(this->tv()[i], i);
	}

	void printChebyshev(const VectorWithOffsetType& phi, SizeType whatTarget) const
	{
		for (SizeType ii = 0; ii < phi.sectors(); ii++) {
			SizeType i = phi.sector(ii);
			printChebyshev(phi, whatTarget, i);
		}
	}

	void printChebyshev(const VectorWithOffsetType& phi, SizeType whatTarget, SizeType i0) const
	{
		SizeType p = this->lrs().super().findPartitionNumber(phi.offset(i0));
		typename ModelType::HamiltonianConnectionType hc(p,
		                                                 BaseType::lrs(),
		                                                 BaseType::model().geometry(),
		                                                 BaseType::model().modelLinks(),
		                                                 this->common().aoe().currentTime(),
		                                                 0);
		typename LanczosSolverType::MatrixType        lanczosHelper(BaseType::model(), hc);

		SizeType         total = phi.effectiveSize(i0);
		TargetVectorType phi2(total);
		phi.extract(phi2, i0);
		TargetVectorType x(total);
		lanczosHelper.matrixVectorProduct(x, phi2);
		PsimagLite::OstringStream                     msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg << "Hamiltonian average at Che-time=" << this->common().aoe().currentTime();
		msg << " for target=" << whatTarget;
		ComplexOrRealType numerator = phi2 * x;
		ComplexOrRealType den       = phi2 * phi2;
		ComplexOrRealType division  = (PsimagLite::norm(den) < 1e-10) ? 0 : numerator / den;
		msg << " sector=" << i0 << " <phi(t)|H|phi(t)>=" << numerator;
		msg << " <phi(t)|phi(t)>=" << den << " " << division;
		progress_.printline(msgg, std::cout);
		tvEnergy_[whatTarget] = PsimagLite::real(division);
	}

	TargetParamsType              tstStruct_;
	const WaveFunctionTransfType& wft_;
	PsimagLite::ProgressIndicator progress_;
	VectorRealType                weight_;
	mutable VectorRealType        tvEnergy_;
	RealType                      gsWeight_;
}; // class TargetingChebyshev
} // namespace Dmrg

#endif
