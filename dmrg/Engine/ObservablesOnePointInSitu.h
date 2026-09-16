// SPDX-FileCopyrightText: Copyright (c) 2009-2014-2021, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 6.]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file ObservablesOnePointInSitu.h
 *
 *  DOC NEEDED FIXME
 */
#ifndef DMRG_OBS_ONE_POINT_IN_SITU_H
#define DMRG_OBS_ONE_POINT_IN_SITU_H

#include "ProgramGlobals.h"

namespace Dmrg {

template <typename VectorType> class ObservablesOnePointInSitu {

public:

	template <typename TargetingType>
	void init(const TargetingType& psi,
	          SizeType             sites) // FIXME : Use ctor instead
	{
		inSitu_.resize(sites, 0);
		for (SizeType i = 0; i < sites; ++i)
			inSitu_[i] = psi.inSitu(i);
	}

	const typename VectorType::value_type& operator()(SizeType i) const { return inSitu_[i]; }

private:

	VectorType inSitu_;
}; // ObservablesOnePointInSitu

} // namespace Dmrg

/*@}*/
#endif // DMRG_OBS_ONE_POINT_IN_SITU_H
