// SPDX-FileCopyrightText: Copyright (c) 2009-2012-2019, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file ParametersGraphene.h
 *
 *  Contains the parameters for the FeAs model and function to read them from a JSON file
 *
 */
#ifndef DMRGPP_PARAMETERS_GRAPHENE_H
#define DMRGPP_PARAMETERS_GRAPHENE_H
#include "ParametersModelBase.h"
#include <PsimagLite/Matrix.h>
#include <stdexcept>
#include <vector>

namespace Dmrg {

//! FeAs Model Parameters
template <typename ComplexOrRealType, typename QnType>
struct ParametersGraphene : public ParametersModelBase<ComplexOrRealType, QnType> {
	// no connections here please!!
	// connections are handled by the geometry

	using RealType = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using BaseType = ParametersModelBase<ComplexOrRealType, QnType>;

	template <typename IoInputType>
	ParametersGraphene(IoInputType& io)
	    : BaseType(io, false)
	{
		io.readline(orbitals, "Orbitals=");
		io.readline(hubbardU, "HubbardU=");
		io.readline(pairHopping, "PairHopping=");
	}

	template <typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType&, SizeType, PsimagLite::String = "") const
	{
		return 0;
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		PsimagLite::String label = label1 + "/ParametersGraphene";
		io.createGroup(label);
		BaseType::write(label, io);
		io.write(label + "/orbitals", orbitals);
		io.write(label + "/hubbardU", hubbardU);
		io.write(label + "/pairHopping", pairHopping);
	}

	SizeType orbitals;
	RealType hubbardU;
	RealType pairHopping;
};
} // namespace Dmrg

/*@}*/
#endif
