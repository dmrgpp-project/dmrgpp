// SPDX-FileCopyrightText: Copyright (c) 2009-2014, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [Lanczos++, Version 1.0.0]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup LanczosPlusPlus */
/*@{*/

/*! \file ParametersFermionSpinless.h
 *
 *  Contains the parameters for the Hubbard model and
 *  function to read them from a JSON file
 *
 */
#ifndef LANCZOSPP_PARAMS_FERMIONSPINLESS_H
#define LANCZOSPP_PARAMS_FERMIONSPINLESS_H
#include <PsimagLite/Vector.h>
#include <stdexcept>

namespace LanczosPlusPlus {
template <typename Field, typename InputType> struct ParametersFermionSpinless {

	ParametersFermionSpinless(InputType& io)
	{
		io.readline(model, "Model=");
		io.read(potentialV, "potentialV");
	}

	ParametersFermionSpinless(const ParametersFermionSpinless&) = delete;

	ParametersFermionSpinless& operator=(const ParametersFermionSpinless&) = delete;

	PsimagLite::String model;
	// Do not include here connection parameters
	// those are handled by the Geometry
	// Onsite potential values, one for each site
	typename PsimagLite::Vector<Field>::Type potentialV;
};

//! Function that prints model parameters to stream os
template <typename FieldType, typename InputType>
std::ostream& operator<<(std::ostream&                                          os,
                         const ParametersFermionSpinless<FieldType, InputType>& params)
{
	PsimagLite::vectorPrint(params.potentialV, "potentialV", os);
	return os;
}
} // namespace LanczosPlusPlus

/*@}*/
#endif
