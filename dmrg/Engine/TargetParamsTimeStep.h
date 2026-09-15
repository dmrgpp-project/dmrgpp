// SPDX-FileCopyrightText: Copyright (c) 2009-2015, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file TargetParamsTimeStep.h
 *
 *  This is a structure to represent the parameters of the TimeStep Evolution
 *  algorithm. Don't add functions to this class because
 *  this class's data is all public
 */
#ifndef TARGET_PARAMS_TIMESTEP_H
#define TARGET_PARAMS_TIMESTEP_H

#include "TargetParamsCommon.h"
#include "TargetParamsTimeVectors.h"
#include <PsimagLite/GetBraOrKet.h>

namespace Dmrg {
// Coordinates reading of TargetSTructure from input file
template <typename ModelType>
class TargetParamsTimeStep : public TargetParamsTimeVectors<ModelType> {

public:

	using TimeVectorParamsType   = TargetParamsTimeVectors<ModelType>;
	using TargetParamsCommonType = TargetParamsCommon<ModelType>;
	using RealType               = typename ModelType::RealType;

	template <typename IoInputter>
	TargetParamsTimeStep(IoInputter& io, PsimagLite::String targeting, const ModelType& model)
	    : TimeVectorParamsType(io, targeting, model)
	    , maxTime_(0)
	{
		try {
			io.readline(maxTime_, "TSPMaxTime=");
		} catch (std::exception&) { }
	}

	RealType maxTime() const override { return maxTime_; }

private:

	RealType maxTime_;
}; // class TargetParamsTimeStep

template <typename ModelType>
inline std::ostream& operator<<(std::ostream& os, const TargetParamsTimeStep<ModelType>& t)
{
	os << "TargetParams.type=TimeStep\n";

	const typename TargetParamsTimeStep<ModelType>::TimeVectorParamsType& tp1 = t;
	os << tp1;

	const typename TargetParamsTimeStep<ModelType>::TargetParamsCommonType& tp = t;
	os << tp;

	if (t.maxTime() > 0)
		os << "TSPMaxTime=" << t.maxTime() << "\n";

	return os;
}
} // namespace Dmrg

/*@}*/
#endif // TARGET_PARAMS_TIMESTEP_H
