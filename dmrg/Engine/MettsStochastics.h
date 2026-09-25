// SPDX-FileCopyrightText: Copyright (c) 2009-2012, 2013, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file MettsStochastics.h
 *
 *  Stochastics (random choices) needed for the METTS algorithm
 *
 */

#ifndef METTS_STOCHASTICS_H
#define METTS_STOCHASTICS_H
#include "Utils.h"
#include <PsimagLite/ProgressIndicator.h>
#include <PsimagLite/TypeToString.h>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace Dmrg {
template <typename ModelType_, typename RngType_> class MettsStochastics {

public:

	using PairType           = std::pair<SizeType, SizeType>;
	using ModelType          = ModelType_;
	using QnType             = typename ModelType::QnType;
	using VectorQnType       = typename QnType::VectorQnType;
	using RealType           = typename ModelType::RealType;
	using LeftRightSuperType = typename ModelType::LeftRightSuperType;
	using HilbertBasisType   = typename ModelType::HilbertBasisType;
	using RngType            = RngType_;
	using LongType           = typename RngType::LongType;
	using VectorRealType     = typename PsimagLite::Vector<RealType>::Type;
	using VectorSizeType     = typename PsimagLite::Vector<SizeType>::Type;

	MettsStochastics(const ModelType& model, int long seed, const VectorSizeType& pure)
	    : model_(model)
	    , rng_(seed)
	    , pure_(pure)
	    , progress_("MettsStochastics")
	    , addedSites_(0)
	{ }

	const ModelType& model() const { return model_; }

	SizeType chooseRandomState(SizeType site) const
	{
		if (site < pure_.size())
			return pure_[site];

		return SizeType(rng_() * model_.hilbertSize(site));
	}

	SizeType chooseRandomState(const VectorRealType& probs) const
	{
		RealType r  = rng_();
		RealType s1 = 0;
		RealType s2 = 0;
		for (SizeType i = 0; i < probs.size(); ++i) {
			s2 = s1 + probs[i];
			if (s1 < r && r <= s2)
				return i;
			s1 = s2;
		}

		PsimagLite::String s(__FILE__);
		s += PsimagLite::String(" ") + ttos(__LINE__) + " " + __FUNCTION__
		    + " Probabilities don't amount to 1\n";
		throw PsimagLite::RuntimeError(s.c_str());
	}

	// call only from INFINITE
	void update(const QnType&                                      qn,
	            const typename PsimagLite::Vector<SizeType>::Type& block1,
	            const typename PsimagLite::Vector<SizeType>::Type& block2,
	            SizeType                                           seed)
	{
		if (addedSites_.size() == 0) {
			pureStates_.resize(block2[block2.size() - 1] + block2.size() + 1);
			initialSetOfPures(seed);
			for (SizeType i = 0; i < block1.size(); i++)
				for (SizeType j = 0; j < block1.size(); j++)
					addedSites_.push_back(block1[i] + j - block1.size());
			for (SizeType i = 0; i < block2.size(); i++)
				for (SizeType j = 0; j < block2.size(); j++)
					addedSites_.push_back(block2[i] + j + block2.size());
		}

		for (SizeType i = 0; i < block1.size(); i++)
			addedSites_.push_back(block1[i]);
		for (SizeType i = 0; i < block2.size(); i++)
			addedSites_.push_back(block2[i]);

		qnVsSize_.resize(addedSites_.size() + 1, QnType::zero());
		qnVsSize_[addedSites_.size()] = qn;
	}

	void setCollapseBasis(typename PsimagLite::Vector<RealType>::Type& collapseBasisWeights,
	                      SizeType                                     site) const
	{
		SizeType nk = model_.hilbertSize(site);
		for (SizeType alpha = 0; alpha < nk; alpha++) {
			RealType randomNumber       = rng_();
			collapseBasisWeights[alpha] = randomNumber;
		}
		RealType norm1 = 1.0 / PsimagLite::norm(collapseBasisWeights);
		collapseBasisWeights *= norm1;
	}

private:

	void initialSetOfPures(LongType)
	{
		for (SizeType i = 0; i < pureStates_.size(); i++)
			pureStates_[i] = SizeType(rng_() * model_.hilbertSize(i));
	}

	const ModelType&                            model_;
	mutable RngType                             rng_;
	const VectorSizeType&                       pure_;
	PsimagLite::ProgressIndicator               progress_;
	typename PsimagLite::Vector<SizeType>::Type pureStates_;
	typename PsimagLite::Vector<SizeType>::Type addedSites_;
	VectorQnType                                qnVsSize_;
}; // class MettsStochastics
} // namespace Dmrg
/*@}*/
#endif // METTS_STOCHASTICS_H
