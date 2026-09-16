// SPDX-FileCopyrightText: Copyright (c) 2009-2012-2021, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file ParametersSu3.h
 *
 *  Contains the parameters for the Heisenberg model and function
 *  to read them from a file
 *
 */
#ifndef DMRG_PARAMS_SU3_H
#define DMRG_PARAMS_SU3_H
#include "../../Engine/ParametersModelBase.h"
#include <PsimagLite/Vector.h>

namespace Dmrg {

template <typename RealType, typename QnType>
struct ParametersSu3 : public ParametersModelBase<RealType, QnType> {

	using BaseType       = ParametersModelBase<RealType, QnType>;
	using VectorRealType = typename PsimagLite::Vector<RealType>::Type;
	// no connectors here, connectors are handled by the geometry
	template <typename IoInputType>
	ParametersSu3(IoInputType& io)
	    : BaseType(io, false)
	{
		io.readline(mass, "Mass=");
		io.readline(p, "Su3RepresentationP=");
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		PsimagLite::String label = label1 + "/ParametersSu3";
		io.createGroup(label);
		BaseType::write(label, io);
		io.write(label + "/mass", mass);
		io.write(label + "/Su3RepresentationP", p);
	}

	RealType mass;
	SizeType p;
};
} // namespace Dmrg

/*@}*/
#endif
