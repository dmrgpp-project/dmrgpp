// SPDX-FileCopyrightText: Copyright (c) 2009-2015, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [by G.A., Oak Ridge National Laboratory]


/** \ingroup DMRG */
/*@{*/

/*! \file ParametersHubbardAncillaExtended.h
 *
 *
 */
#ifndef ParametersHubbardAncillaExtended_H
#define ParametersHubbardAncillaExtended_H
#include "../FeAsModel/ParametersModelFeAs.h"
#include "ParametersModelBase.h"

namespace Dmrg {
//! Hubbard Model Parameters
template <typename RealType, typename QnType>
struct ParametersHubbardAncillaExtended : ParametersModelFeAs<RealType, QnType> {

	using BaseType = ParametersModelFeAs<RealType, QnType>;

	template <typename IoInputType>
	ParametersHubbardAncillaExtended(IoInputType& io)
	    : BaseType(io)
	{ }

	template <typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType&, SizeType, PsimagLite::String = "") const
	{
		return 0;
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		PsimagLite::String label = label1 + "/ParametersHubbardAncillaExtended";
		io.createGroup(label);
		BaseType::write(label, io);
	}
};

} // namespace Dmrg

/*@}*/
#endif
