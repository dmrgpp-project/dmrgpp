// SPDX-FileCopyrightText: Copyright (c) 2009-2014, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file TargetHelper.h
 *
 *  FIXME
 */
#ifndef TARGET_HELPER_H
#define TARGET_HELPER_H

namespace Dmrg {

template <typename ModelType_, typename WaveFunctionTransfType_> class TargetHelper {

public:

	using ModelType              = ModelType_;
	using WaveFunctionTransfType = WaveFunctionTransfType_;
	using RealType               = typename ModelType::RealType;
	using ModelHelperType        = typename ModelType::ModelHelperType;
	using LeftRightSuperType     = typename ModelHelperType::LeftRightSuperType;
	using InputValidatorType     = typename ModelType::InputValidatorType;

	TargetHelper(const LeftRightSuperType&     lrs,
	             const ModelType&              model,
	             const WaveFunctionTransfType& wft)
	    : lrs_(lrs)
	    , model_(model)
	    , wft_(wft)
	{ }

	const LeftRightSuperType& lrs() const { return lrs_; }

	const ModelType& model() const { return model_; }

	const WaveFunctionTransfType& wft() const { return wft_; }

	bool withLegacyBugs() const { return model_.params().options.isSet("keepLegacyBugs"); }

private:

	const LeftRightSuperType&     lrs_;
	const ModelType&              model_;
	const WaveFunctionTransfType& wft_;
}; // TargetHelper

} // namespace Dmrg

/*@}*/
#endif // TARGET_HELPER_H
