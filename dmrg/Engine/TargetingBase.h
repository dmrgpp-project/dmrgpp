// SPDX-FileCopyrightText: Copyright (c) 2014, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file TargetingBase.h
 *
 *
 *
 */

#ifndef TARGETING_BASE_H
#define TARGETING_BASE_H
#include "Checkpoint.h"
#include "Intent.h"
#include "OneSiteSpaces.hh"
#include "TargetHelper.h"
#include "TargetParamsBase.h"
#include "TargetingCommon.h"
#include "Wft/WaveFunctionTransfFactory.h"
#include <PsimagLite/Io/IoSelector.h>
#include <iostream>

namespace Dmrg {

template <typename LanczosSolverType_, typename VectorWithOffsetType_> class TargetingBase {

public:

	using LanczosSolverType       = LanczosSolverType_;
	using VectorWithOffsetType    = VectorWithOffsetType_;
	using MatrixVectorType        = typename LanczosSolverType::MatrixType;
	using ModelType               = typename MatrixVectorType::ModelType;
	using RealType                = typename ModelType::RealType;
	using ParametersType          = typename ModelType::ParametersType;
	using OptionsType             = typename ParametersType::OptionsType;
	using ParametersForSolverType = PsimagLite::ParametersForSolver<RealType>;
	using ModelHelperType         = typename ModelType::ModelHelperType;
	using LeftRightSuperType      = typename ModelHelperType::LeftRightSuperType;
	using BasisWithOperatorsType  = typename LeftRightSuperType::BasisWithOperatorsType;
	using SparseMatrixType        = typename BasisWithOperatorsType::SparseMatrixType;
	using OperatorType            = typename BasisWithOperatorsType::OperatorType;
	using BasisType               = typename BasisWithOperatorsType::BasisType;
	using ComplexOrRealType       = typename SparseMatrixType::value_type;
	using BlockType               = typename BasisType::BlockType;
	using QnType                  = typename BasisType::QnType;
	using OneSiteSpacesType       = OneSiteSpaces<ModelType>;
	using WaveFunctionTransfType  = WaveFunctionTransfFactory<LeftRightSuperType,
	                                                          VectorWithOffsetType,
	                                                          OptionsType,
	                                                          OneSiteSpacesType>;
	using VectorType              = typename VectorWithOffsetType::VectorType;
	using TargetVectorType        = VectorType;
	using TargetParamsType        = TargetParamsBase<ModelType>;
	using TargetHelperType        = TargetHelper<ModelType, WaveFunctionTransfType>;
	using TargetingCommonType
	    = TargetingCommon<TargetHelperType, VectorWithOffsetType, LanczosSolverType>;
	typedef
	    typename TargetingCommonType::ApplyOperatorExpressionType ApplyOperatorExpressionType;
	using VectorOperatorType     = typename PsimagLite::Vector<OperatorType>::Type;
	using StageEnumType          = typename ApplyOperatorExpressionType::StageEnumType;
	using DmrgSerializerType     = typename ApplyOperatorExpressionType::DmrgSerializerType;
	using VectorSizeType         = typename PsimagLite::Vector<SizeType>::Type;
	using VectorRealType         = typename PsimagLite::Vector<RealType>::Type;
	using VectorVectorType       = typename PsimagLite::Vector<TargetVectorType>::Type;
	using VectorVectorVectorType = typename PsimagLite::Vector<VectorVectorType>::Type;
	using VectorVectorVectorWithOffsetType =
	    typename TargetingCommonType::VectorVectorVectorWithOffsetType;
	using CheckpointType = Checkpoint<ModelType, WaveFunctionTransfType>;

	TargetingBase(const LeftRightSuperType&     lrs,
	              const CheckpointType&         checkPoint,
	              const WaveFunctionTransfType& wft,
	              SizeType                      indexNoAdvance)
	    : lrs_(lrs)
	    , model_(checkPoint.model())
	    , commonTargeting_(lrs, checkPoint, wft, indexNoAdvance)
	{
		Intent<ModelType> intent(model_);
		intent.check();

		const SizeType nexcited = model_.params().numberOfExcited;

		if (nexcited == 1)
			return; // EARLY EXIT

		PsimagLite::String msg = "nexcited = " + ttos(nexcited) + " > 1 is experimental\n";
		std::cerr << msg;
		std::cout << msg;
	}

	TargetingBase(const TargetingBase&) = delete;

	TargetingBase& operator=(const TargetingBase&) = delete;

	virtual ~TargetingBase() { }

	virtual void postCtor() { commonTargeting_.postCtor(sites(), targets()); }

	virtual SizeType sites() const = 0;

	virtual SizeType targets() const = 0;

	virtual RealType gsWeight() const = 0;

	virtual RealType weight(SizeType i) const = 0;

	virtual void evolve(const VectorRealType&         energies,
	                    ProgramGlobals::DirectionEnum direction,
	                    const BlockType&              block1,
	                    const BlockType&              block2,
	                    SizeType                      loopNumber)
	    = 0;

	virtual void read(typename TargetingCommonType::IoInputType&, PsimagLite::String) = 0;

	virtual void
	write(const VectorSizeType&, PsimagLite::IoSelector::Out&, PsimagLite::String) const
	    = 0;

	// virtuals with default implementation

	virtual bool includeGroundStage() const { return true; }

	virtual void
	set(VectorVectorVectorType& inV, const VectorSizeType& sectors, const BasisType& someBasis)
	{
		commonTargeting_.aoeNonConst().setPsi(
		    inV, sectors, someBasis, model_.params().numberOfExcited);
	}

	virtual void updateOnSiteForCorners(BasisWithOperatorsType& basisWithOps) const
	{
		if (BasisWithOperatorsType::useSu2Symmetry())
			return;

		BlockType X = basisWithOps.block();

		if (X.size() != 1)
			return;

		if (X[0] != 0 && X[0] != lrs_.super().block().size() - 1)
			return;

		basisWithOps.setOneSite(X, model_, commonTargeting_.time());
	}

	virtual bool end() const { return false; }

	virtual SizeType size() const
	{
		if (commonTargeting_.aoe().allStages(StageEnumType::DISABLED))
			return 0;
		return commonTargeting_.aoe().tvs();
	}

	virtual RealType normSquared(SizeType i) const { return commonTargeting_.normSquared(i); }

	virtual void initPsi(SizeType nsectors, SizeType nexcited)
	{
		commonTargeting_.aoeNonConst().initPsi(nsectors, nexcited);
	}

	// legacy thing for vectorwithoffsets
	virtual void initialGuess(VectorVectorType&        initialVector,
	                          const OneSiteSpacesType& oneSiteSpaces,
	                          bool                     noguess,
	                          const VectorSizeType&    compactedWeights,
	                          const VectorSizeType&    sectors,
	                          const BasisType&         basis) const
	{
		if (!model_.params().options.isSet("vectorwithoffsets"))
			err("FATAL: Wrong execution path\n");

		const VectorWithOffsetType& psi00
		    = commonTargeting_.aoe().ensureOnlyOnePsi("initialGuess");
		VectorWithOffsetType vwo(compactedWeights, sectors, basis);
		commonTargeting_.initialGuess(vwo, psi00, oneSiteSpaces, noguess);
		const SizeType n = vwo.sectors();
		initialVector.resize(n);
		for (SizeType i = 0; i < n; ++i)
			vwo.extract(initialVector[i], vwo.sector(i));
	}

	virtual void initialGuess(VectorType&              initialVector,
	                          const OneSiteSpacesType& oneSiteSpaces,
	                          bool                     noguess,
	                          const VectorSizeType&    compactedWeights,
	                          const VectorSizeType&    sectors,
	                          SizeType                 sectorIndex,
	                          SizeType                 excited,
	                          const BasisType&         basis) const
	{

		if (model_.params().options.isSet("vectorwithoffsets"))
			err("FATAL: Wrong execution path\n");

		const VectorVectorVectorWithOffsetType& psi = commonTargeting_.aoe().psiConst();
		const SizeType                          nsectors = psi.size();

		if (nsectors != compactedWeights.size())
			err("initialGuess compactedWeights\n");

		const SizeType numberOfExcited = psi[sectorIndex].size();

		if (excited > numberOfExcited)
			err("initialGuess, excited=" + ttos(excited) + " > " + ttos(numberOfExcited)
			    + "\n");

		SizeType start = 0;
		SizeType end   = numberOfExcited;
		if (excited < numberOfExcited) {
			start = excited;
			end   = excited + 1;
		}

		assert(sectorIndex < psi.size());
		assert(sectorIndex < compactedWeights.size());
		assert(sectorIndex < sectors.size());
		for (SizeType e = start; e < end; ++e) {

			VectorWithOffsetType vwo(
			    compactedWeights[sectorIndex], sectors[sectorIndex], basis);

			assert(e < psi[sectorIndex].size());
			if (psi[sectorIndex][e] == nullptr) {
				noguess = true;
			}

			commonTargeting_.initialGuess(
			    vwo, *(psi[sectorIndex][e]), oneSiteSpaces, noguess);

			VectorType tmpVector;
			vwo.extract(tmpVector, vwo.sector(0));
			if (e == start)
				initialVector = tmpVector;
			else
				initialVector += tmpVector;
		}
	}

	// non-virtual below

	const ModelType& model() const { return model_; }

	const VectorVectorVectorWithOffsetType& psiConst() const
	{
		return commonTargeting_.aoe().psiConst();
	}

	const VectorWithOffsetType& operator()(SizeType i) const
	{
		return commonTargeting_.aoe().targetVectors(i);
	}

	RealType time() const { return commonTargeting_.time(); }

	const typename VectorWithOffsetType::value_type& inSitu(SizeType i) const
	{
		return commonTargeting_.inSitu(i);
	}

	const LeftRightSuperType& lrs() const { return lrs_; }

	static PsimagLite::String buildPrefix(PsimagLite::IoSelector::Out& io, SizeType counter)
	{
		PsimagLite::String prefix("TargetingCommon");
		using SerializerType = PsimagLite::IoSelector::Out::Serializer;
		if (counter == 0)
			io.createGroup(prefix);

		io.write(counter + 1,
		         prefix + "/Size",
		         (counter == 0) ? SerializerType::NO_OVERWRITE
		                        : SerializerType::ALLOW_OVERWRITE);

		prefix += ("/" + ttos(counter));

		io.createGroup(prefix);
		return prefix;
	}

	void multiSitePush(DmrgSerializerType const* ds) const
	{
		commonTargeting_.aoe().multiSitePush(ds);
	}

protected:

	TargetingCommonType& common() { return commonTargeting_; }

	const TargetingCommonType& common() const { return commonTargeting_; }

	VectorWithOffsetType& tvNonConst(SizeType ind)
	{
		return commonTargeting_.aoeNonConst().targetVectorsNonConst(ind);
	}

	const VectorWithOffsetType& tv(SizeType ind) const
	{
		return commonTargeting_.aoe().targetVectors(ind);
	}

	SizeType numberOfTvs() const { return commonTargeting_.aoe().tvs(); }

private:

	const LeftRightSuperType& lrs_;
	const ModelType&          model_;
	TargetingCommonType       commonTargeting_;
}; // class TargetingBase

} // namespace Dmrg
/*@}*/
#endif
