// SPDX-FileCopyrightText: Copyright (c) 2009-2014, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [by G.A., Oak Ridge National Laboratory]


#ifndef TARGETING_TIMESTEP_H
#define TARGETING_TIMESTEP_H

#include "BlockDiagonalMatrix.h"
#include "ProgramGlobals.h"
#include "TargetParamsTimeStep.h"
#include "TargetingBase.h"
#include "TimeVectorsKrylov.h"
#include "TimeVectorsRungeKutta.h"
#include "TimeVectorsSuzukiTrotter.h"
#include <PsimagLite/ParametersForSolver.h>
#include <PsimagLite/PredicateAwesome.h>
#include <PsimagLite/ProgressIndicator.h>
#include <iostream>

namespace Dmrg {

template <typename LanczosSolverType_, typename VectorWithOffsetType_>
class TargetingTimeStep : public TargetingBase<LanczosSolverType_, VectorWithOffsetType_> {

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
	using CheckpointType         = typename BaseType::CheckpointType;
	using PairType               = std::pair<SizeType, SizeType>;
	using OptionsType            = typename BaseType::OptionsType;
	using MatrixVectorType       = typename BaseType::MatrixVectorType;
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
	using StageEnumType          = typename TargetingCommonType::StageEnumType;

	TargetingTimeStep(const LeftRightSuperType&     lrs,
	                  const CheckpointType&         checkPoint,
	                  const WaveFunctionTransfType& wft,
	                  const QnType&,
	                  InputValidatorType& ioIn,
	                  PsimagLite::String  targeting)
	    : BaseType(lrs, checkPoint, wft, 0)
	    , tstStruct_(ioIn, targeting, checkPoint.model())
	    , wft_(wft)
	    , progress_(targeting)
	    , weight_(tstStruct_.times().size())
	    , tvEnergy_(tstStruct_.times().size(), 0.0)
	    , gsWeight_(tstStruct_.gsWeight())
	{
		if (!wft.isEnabled())
			err("TST needs an enabled wft\n");
		if (tstStruct_.sites() == 0)
			err("TST needs at least one TSPSite\n");

		RealType tau = tstStruct_.tau();
		RealType sum = 0;
		SizeType n   = tstStruct_.times().size();

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

	bool end() const override
	{
		return (tstStruct_.maxTime() != 0
		        && this->common().aoe().timeVectors().time() >= tstStruct_.maxTime());
	}

	void read(typename TargetingCommonType::IoInputType& io, PsimagLite::String prefix) override
	{
		this->common().readGSandNGSTs(io, prefix, "TimeStep");
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

		SizeType numberOfSites    = this->lrs().super().block().size();
		bool     doBorderIfBorder = (site < 1 || site >= numberOfSites - 1);

		if (doBorderIfBorder) {
			if (loopNumber >= this->model().params().finiteLoop.size() - 1) {
				if (direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) {
					if (site >= numberOfSites - 1) {
						this->common().cocoon(block1, direction, false);
						return;
					}
				} else {
					if (site < 1) {
						this->common().cocoon(block1, direction, false);
						return;
					}
				}
			}
			this->common().cocoon(block1, direction, false);
		}

		this->common().aoeNonConst().getPhi(
		    &phiNew, Eg, direction, site, loopNumber, tstStruct_);

		PairType startEnd(0, tstStruct_.times().size());
		bool     allOperatorsApplied
		    = (this->common().aoe().noStageIs(StageEnumType::DISABLED)
		       && this->common().aoe().noStageIs(StageEnumType::OPERATOR));

		VectorSizeType indices(startEnd.second - startEnd.first);
		for (SizeType i = 0; i < indices.size(); ++i)
			indices[i] = i + startEnd.first;

		static const bool isLastCall = true;
		this->common().aoeNonConst().calcTimeVectors(
		    indices,
		    Eg,
		    phiNew,
		    direction,
		    allOperatorsApplied,
		    false, // don't wft or advance indices[0]
		    block1,
		    isLastCall);

		this->common().cocoon(block1, direction, false);

		PsimagLite::String predicate = this->model().params().printHamiltonianAverage;
		const SizeType     center    = this->model().superGeometry().numberOfSites() / 2;
		PsimagLite::replaceAll(predicate, "c", ttos(center));
		PsimagLite::PredicateAwesome<> pAwesome(predicate);
		assert(block1.size() > 0);
		if (pAwesome.isTrue("s", block1[0]))
			printEnergies(); // in-situ

		const OptionsType& options = this->model().params().options;
		bool               normalizeTimeVectors
		    = (options.isSet("normalizeTimeVectors") || options.isSet("TargetingAncilla"));

		if (options.isSet("neverNormalizeVectors"))
			normalizeTimeVectors = false;

		if (normalizeTimeVectors)
			this->common().normalizeTimeVectors();

		this->common().printNormsAndWeights(gsWeight_, weight_);
	}

	void printEnergies() const
	{
		for (SizeType i = 0; i < this->common().aoe().tvs(); i++)
			printEnergies(this->tv(i), i);
	}

	void printEnergies(const VectorWithOffsetType& phi, SizeType whatTarget) const
	{
		for (SizeType ii = 0; ii < phi.sectors(); ii++) {
			SizeType i = phi.sector(ii);
			printEnergies(phi, whatTarget, i);
		}
	}

	void printEnergies(const VectorWithOffsetType& phi, SizeType whatTarget, SizeType i0) const
	{
		const SizeType p = this->lrs().super().findPartitionNumber(phi.offset(i0));
		typename ModelHelperType::Aux                 aux(p, BaseType::lrs());
		typename ModelType::HamiltonianConnectionType hc(
		    BaseType::lrs(),
		    ModelType::modelLinks(),
		    this->common().aoe().timeVectors().time(),
		    BaseType::model().superOpHelper(),
		    BaseType::model().ioIn());
		typename LanczosSolverType::MatrixType lanczosHelper(BaseType::model(), hc, aux);

		const SizeType   total = phi.effectiveSize(i0);
		TargetVectorType phi2(total);
		phi.extract(phi2, i0);
		TargetVectorType x(total);
		lanczosHelper.matrixVectorProduct(x, phi2);
		PsimagLite::OstringStream                     msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg << "Hamiltonian average at time=" << this->common().aoe().timeVectors().time();
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
}; // class TargetingTimeStep
} // namespace Dmrg

#endif
