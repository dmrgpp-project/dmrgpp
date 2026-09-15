// SPDX-FileCopyrightText: Copyright (c) 2009-2013, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file TargetParamsGroundState.h
 *
 * There's nothing to see here.
 */
#ifndef TARGET_PARAMS_GROUNDSTATE_H
#define TARGET_PARAMS_GROUNDSTATE_H
#include "TargetParamsBase.h"

namespace Dmrg {
// Coordinates reading of TargetSTructure from input file
template <typename ModelType> class TargetParamsGroundState : public TargetParamsBase<ModelType> {

	using BaseType = TargetParamsBase<ModelType>;

public:

	using RealType = typename ModelType::RealType;

	TargetParamsGroundState(PsimagLite::String targeting)
	    : BaseType(targeting)
	{ }

	SizeType sites() const override { return 0; }

	SizeType sectorIndex() const override
	{
		throw PsimagLite::RuntimeError("sectorIndex called for gs\n");
	}

	SizeType levelIndex() const override
	{
		throw PsimagLite::RuntimeError("levelIndex called for gs\n");
	}

	void write(PsimagLite::String, PsimagLite::IoSerializer&) const
	{
		// this function intentionally left blank
	}
}; // class TargetParamsGroundState

template <typename ModelType>
inline std::ostream& operator<<(std::ostream& os, const TargetParamsGroundState<ModelType>&)
{
	os << "TargetParams.type=GroundState\n";
	return os;
}
} // namespace Dmrg

/*@}*/
#endif // TARGET_PARAMS_GROUNDSTATE_H
