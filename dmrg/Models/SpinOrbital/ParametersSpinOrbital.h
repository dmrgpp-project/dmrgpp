// SPDX-FileCopyrightText: Copyright (c) 2009-2021, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [by G.A., Oak Ridge National Laboratory]


/** \ingroup DMRG */
/*@{*/

/*! \file ParametersSpinOrbital.h
 *
 *  Contains the parameters for the Heisenberg model and function
 *  to read them from a file
 *
 */
#ifndef DMRG_PARAMS_SPIN_ORBITAL_H
#define DMRG_PARAMS_SPIN_ORBITAL_H
#include "../../Engine/ParametersModelBase.h"
#include <PsimagLite/Vector.h>

namespace Dmrg {
template <typename RealType, typename QnType>
struct ParametersSpinOrbital : public ParametersModelBase<RealType, QnType> {

	using BaseType       = ParametersModelBase<RealType, QnType>;
	using VectorRealType = typename PsimagLite::Vector<RealType>::Type;
	// no connectors here, connectors are handled by the geometry
	template <typename IoInputType>
	ParametersSpinOrbital(IoInputType& io)
	    : BaseType(io, false)
	{
		io.readline(twiceS, "SpinTwiceS=");
		io.readline(twiceL, "OrbitalTwiceS=");

		io.readline(lambda1, "LambdaOne=");
		io.readline(lambda2, "LambdaTwo=");
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		PsimagLite::String label = label1 + "/ParametersSpinOrbital";
		io.createGroup(label);
		BaseType::write(label, io);
		io.write(label + "/SpinTwiceS", twiceS);
		io.write(label + "/OrbitalTwiceS", twiceL);
		io.write(label + "/LambdaOne", lambda1);
		io.write(label + "/LambdaTwo", lambda2);
	}

	//! Function that prints model parameters to stream os
	friend std::ostream& operator<<(std::ostream& os, const ParametersSpinOrbital& parameters)
	{
		os << "SpinTwiceS=" << parameters.twiceS << "\n";
		os << "OrbitalTwiceS=" << parameters.twiceL << "\n";
		os << "LambdaOne=" << parameters.lambda1 << "\n";
		os << "LambdaTwo=" << parameters.lambda2 << "\n";
		os << parameters.targetQuantum;
		return os;
	}

	SizeType twiceS;
	SizeType twiceL;
	RealType lambda1;
	RealType lambda2;
};
} // namespace Dmrg

/*@}*/
#endif
