// SPDX-FileCopyrightText: Copyright (c) 2009-2015, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file ParametersTjAncillaC.h
 *
 *  Contains the parameters for the Hubbard model and function to
 *  read them from a file
 *
 */
#ifndef DMRG_PARAMS_TJ_ANCILLAC_H
#define DMRG_PARAMS_TJ_ANCILLAC_H
#include "ParametersModelBase.h"

namespace Dmrg {
//! TjAncillaC Model Parameters
// no connectors here, connectors are handled by the geometry
template <typename RealType, typename QnType>
struct ParametersTjAncillaC : public ParametersModelBase<RealType, QnType> {

	using BaseType = ParametersModelBase<RealType, QnType>;

	template <typename IoInputType>
	ParametersTjAncillaC(IoInputType& io)
	    : BaseType(io, false)
	{
		io.read(potentialV, "potentialV");
	}

	template <typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType&, SizeType, PsimagLite::String = "") const
	{
		return 0;
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		PsimagLite::String label = label1 + "/ParametersTjAncillaC";
		io.createGroup(label);
		BaseType::write(label, io);
		io.write(label + "/potentialV", potentialV);
	}

	//! Function that prints model parameters to stream os
	friend std::ostream& operator<<(std::ostream& os, const ParametersTjAncillaC& parameters)
	{
		os << "potentialV\n";
		os << parameters.potentialV;
		return os;
	}

	// Do not include here connection parameters
	typename PsimagLite::Vector<RealType>::Type potentialV;
};
} // namespace Dmrg

/*@}*/
#endif
