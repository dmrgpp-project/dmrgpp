// SPDX-FileCopyrightText: Copyright (c) 2009-2015, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [Lanczos++, Version 3.0]
// [by G.A., Oak Ridge National Laboratory]


/** \ingroup LanczosPlusPlus */
/*@{*/

/*! \file ParametersHeisenberg.h
 *
 *  Contains the parameters for the Tj1Orb model
 *
 */
#ifndef LANCZOS_PARAMS_HEISENBERG_H
#define LANCZOS_PARAMS_HEISENBERG_H
#include <PsimagLite/Vector.h>
#include <stdexcept>

namespace LanczosPlusPlus {
//! Heisenberg Model Parameters
template <typename RealType, typename InputType> struct ParametersHeisenberg {

	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

	ParametersHeisenberg(InputType& io)
	{
		io.readline(twiceTheSpin, "HeisenbergTwiceS=");

		try {
			io.read(magneticField, "MagneticField");
		} catch (std::exception&) { }

		try {
			io.read(anisotropy, "AnisotropyD");
		} catch (std::exception&) { }
	}

	// Do not include here connection parameters
	// those are handled by the Geometry
	SizeType       twiceTheSpin;
	VectorRealType magneticField;
	VectorRealType anisotropy;
};

//! Function that prints model parameters to stream os
template <typename RealType, typename InputType>
std::ostream& operator<<(std::ostream&                                    os,
                         const ParametersHeisenberg<RealType, InputType>& parameters)
{
	os << "MagneticField=" << parameters.magneticField << "\n";
	os << "AnisotropyD=" << parameters.anisotropy << "\n";
	os << "HeisenbergTwiceS=" << parameters.twiceTheSpin << "\n";
	return os;
}
} // namespace LanczosPlusPlus

/*@}*/
#endif
