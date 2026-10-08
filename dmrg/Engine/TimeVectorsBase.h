// SPDX-FileCopyrightText: Copyright (c) 2009-2012, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 2.0.0]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file TimeVectorsBase.h
 *
 *
 */

#ifndef TIME_VECTORS_BASE
#define TIME_VECTORS_BASE
#include "ProgramGlobals.h"
#include "Wft/WftHelper.h"
#include <PsimagLite/Vector.h>
#include <iostream>

namespace Dmrg {

template <typename TargetParamsType,
          typename ModelType,
          typename WaveFunctionTransfType,
          typename LanczosSolverType,
          typename VectorWithOffsetType>
class TimeVectorsBase {

public:

	using ModelHelperType            = typename ModelType::ModelHelperType;
	using LeftRightSuperType         = typename ModelHelperType::LeftRightSuperType;
	using RealType                   = typename TargetParamsType::RealType;
	using VectorRealType             = typename PsimagLite::Vector<RealType>::Type;
	using PairType                   = std::pair<SizeType, SizeType>;
	using VectorVectorWithOffsetType = typename PsimagLite::Vector<VectorWithOffsetType>::Type;
	using VectorSizeType             = PsimagLite::Vector<SizeType>::Type;
	using WftHelperType = WftHelper<ModelType, VectorWithOffsetType, WaveFunctionTransfType>;

	TimeVectorsBase(const ModelType&              model,
	                const LeftRightSuperType&     lrs,
	                const WaveFunctionTransfType& wft,
	                PsimagLite::String            name)
	    : wftHelper_(model, lrs, wft)
	    , name_(name)
	    , time_(0.0)
	    , currentTimeStep_(0)
	{ }

	struct ExtraData {

		ExtraData(ProgramGlobals::DirectionEnum dir_,
		          bool                          allOperatorsApplied_,
		          bool                          wftAndAdvanceIfNeeded_,
		          VectorSizeType                block_,
		          bool                          isLastCall_)
		    : dir(dir_)
		    , allOperatorsApplied(allOperatorsApplied_)
		    , wftAndAdvanceIfNeeded(wftAndAdvanceIfNeeded_)
		    , block(block_)
		    , isLastCall(isLastCall_)
		{ }

		ProgramGlobals::DirectionEnum      dir;
		bool                               allOperatorsApplied;
		bool                               wftAndAdvanceIfNeeded;
		PsimagLite::Vector<SizeType>::Type block;
		bool                               isLastCall;
	};

	virtual void calcTimeVectors(const VectorSizeType&,
	                             RealType,
	                             const VectorWithOffsetType&,
	                             const ExtraData&)
	{
		err("calcTimeVectors: unimplemented in this base class\n");
	}

	virtual ~TimeVectorsBase() { }

	virtual void timeHasAdvanced() { }

	RealType time() const { return time_; }

	SizeType currentTimeStep() const { return currentTimeStep_; }

	void advanceCurrentTime(RealType tau) { time_ += tau; }

	void advanceCurrentTimeStep() { ++currentTimeStep_; }

	void setCurrentTimeStep(SizeType t) { currentTimeStep_ = t; }

	void setCurrentTime(RealType t) { time_ = t; }

	bool isBase() const { return (name_ == "base"); }

protected:

	const WftHelperType& wftHelper() const { return wftHelper_; }

private:

	WftHelperType      wftHelper_;
	PsimagLite::String name_;
	RealType           time_;
	SizeType           currentTimeStep_;
}; // class TimeVectorsBase
} // namespace Dmrg
/*@}*/
#endif
