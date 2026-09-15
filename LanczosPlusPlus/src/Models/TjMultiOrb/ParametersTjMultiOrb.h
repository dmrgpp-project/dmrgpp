// SPDX-FileCopyrightText: Copyright (c) 2009-2012, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [Lanczos++, Version 2.0.0]
// [by G.A., Oak Ridge National Laboratory]


/** \ingroup LanczosPlusPlus */
/*@{*/

/*! \file ParametersTjMultiOrb.h
 *
 *  Contains the parameters for the Tj1Orb model
 *
 */
#ifndef LANCZOS_PARAMS_TJ_MULTIORB_H
#define LANCZOS_PARAMS_TJ_MULTIORB_H
#include <PsimagLite/Vector.h>
#include <stdexcept>

namespace LanczosPlusPlus {
//! Tj1Orb Model Parameters
template <typename Field, typename InputType> struct ParametersTjMultiOrb {

	ParametersTjMultiOrb(InputType& io)
	    : reinterpretAndTruncate(0)
	{
		try {
			io.read(potentialV, "potentialV");
		} catch (std::exception&) { }

		try {
			io.readline(reinterpretAndTruncate, "JHundInfinity=");
		} catch (std::exception&) { }

		io.read(orbitals, "Orbitals=");
		if (orbitals != 2 && reinterpretAndTruncate) {
			throw PsimagLite::RuntimeError(
			    "JHundInfinity=1 only possible for orbitals==2\n");
		}
	}

	// Do not include here connection parameters
	// those are handled by the Geometry
	SizeType                                 orbitals;
	typename PsimagLite::Vector<Field>::Type potentialV;
	SizeType                                 reinterpretAndTruncate;
};

//! Function that prints model parameters to stream os
template <typename FieldType, typename InputType>
std::ostream& operator<<(std::ostream& os, const ParametersTjMultiOrb<FieldType, InputType>& p)
{
	os << "Orbitals=" << p.orbitals << "\n";
	os << "JHundInfinity=" << p.reinterpretAndTruncate << "\n";
	return os;
}
} // namespace LanczosPlusPlus

/*@}*/
#endif
