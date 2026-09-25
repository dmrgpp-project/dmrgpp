// SPDX-FileCopyrightText: Copyright (c) 2009-2014, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [Lanczos++, Version 1.0.0]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup LanczosPlusPlus */
/*@{*/

/*! \file ParametersModelHubbard.h
 *
 *  Contains the parameters for the Hubbard model and
 *  function to read them from a JSON file
 *
 */
#ifndef LANCZOS_PARAMS_HUBBARD_H
#define LANCZOS_PARAMS_HUBBARD_H
#include <PsimagLite/Vector.h>
#include <stdexcept>

namespace LanczosPlusPlus {
//! Hubbard Model Parameters
template <typename Field, typename InputType> struct ParametersModelHubbard {

	ParametersModelHubbard(InputType& io)
	{
		io.readline(model, "Model=");
		io.read(hubbardU, "hubbardU");
		io.read(potentialV, "potentialV");
		timeFactor = 0;
		try {
			io.read(potentialT, "PotentialT");
			io.readline(timeFactor, "timeFactor=");
		} catch (std::exception& e) { }
	}

	ParametersModelHubbard(const ParametersModelHubbard&) = delete;

	ParametersModelHubbard& operator=(const ParametersModelHubbard&) = delete;

	PsimagLite::String model;
	// Do not include here connection parameters
	// those are handled by the Geometry
	// Hubbard U values (one for each site)
	typename PsimagLite::Vector<Field>::Type hubbardU;
	// Onsite potential values, one for each site
	typename PsimagLite::Vector<Field>::Type potentialV;
	typename PsimagLite::Vector<Field>::Type potentialT;
	Field                                    timeFactor;
};

//! Function that prints model parameters to stream os
template <typename FieldType, typename InputType>
std::ostream& operator<<(std::ostream&                                       os,
                         const ParametersModelHubbard<FieldType, InputType>& params)
{
	PsimagLite::vectorPrint(params.hubbardU, "hubbardU", os);
	PsimagLite::vectorPrint(params.potentialV, "potentialV", os);
	return os;
}
} // namespace LanczosPlusPlus

/*@}*/
#endif
