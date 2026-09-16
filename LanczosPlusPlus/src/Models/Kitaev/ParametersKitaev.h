// SPDX-FileCopyrightText: Copyright (c) 2009-2015-2018, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [Lanczos++, Version 1.]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup LanczosPlusPlus */
/*@{*/

/*! \file ParametersKitaev.h
 *
 *  Contains the parameters for the Kitaev model
 *
 */
#ifndef LANCZOS_PARAMS_KITAEV_H
#define LANCZOS_PARAMS_KITAEV_H
#include <PsimagLite/Vector.h>
#include <stdexcept>

namespace LanczosPlusPlus {
//! Kitaev Model Parameters
template <typename RealType, typename InputType> struct ParametersKitaev {

	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

	ParametersKitaev(InputType& io)
	{
		try {
			io.read(magneticField, "MagneticField");
		} catch (std::exception&) { }
	}

	// Do not include here connection parameters
	// those are handled by the Geometry
	VectorRealType magneticField;
};

//! Function that prints model parameters to stream os
template <typename RealType, typename InputType>
std::ostream& operator<<(std::ostream& os, const ParametersKitaev<RealType, InputType>& parameters)
{
	os << "MagneticField=" << parameters.magneticField << "\n";
	return os;
}
} // namespace LanczosPlusPlus

/*@}*/
#endif
