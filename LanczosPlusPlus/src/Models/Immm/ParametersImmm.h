// SPDX-FileCopyrightText: Copyright (c) 2009, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 2.0.0]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file ParametersImmm.h
 *
 *  Contains the parameters for the FeAs model and function to read them from a JSON file
 *
 */
#ifndef LANCZOS_PARAMS_IMMM_H
#define LANCZOS_PARAMS_IMMM_H

namespace LanczosPlusPlus {
//! FeAs Model Parameters
template <typename Field> struct ParametersImmm {
	// no connections here please!!
	// connections are handled by the geometry

	template <typename IoInputType> ParametersImmm(IoInputType& io)
	{

		io.read(hubbardU, "hubbardU");
		io.read(potentialV, "potentialV");
		useReflectionSymmetry = 0;
		nOfElectrons          = 0;
		try {
			io.readline(useReflectionSymmetry, "UseReflectionSymmetry=");
		} catch (std::exception& e) { }
	}

	// Hubbard U values (one for each site)
	typename PsimagLite::Vector<Field>::Type hubbardU;
	// Onsite potential values, one for each site
	typename PsimagLite::Vector<Field>::Type potentialV;
	// use reflection symmetry if set to 1
	SizeType useReflectionSymmetry;
	// target number of electrons  in the system
	int nOfElectrons;
	// target density
	// Field density;
}; // struct ParametersImmm

//! Function that prints model parameters to stream os
template <typename FieldType>
std::ostream& operator<<(std::ostream& os, const ParametersImmm<FieldType>& parameters)
{
	// os<<"parameters.density="<<parameters.density<<"\n";
	os << "hubbardU\n";
	os << parameters.hubbardU;
	os << "potentialV\n";
	os << parameters.potentialV;
	os << "UseReflectionSymmetry=" << parameters.useReflectionSymmetry << "\n";
	return os;
}
} // namespace Dmrg

/*@}*/
#endif // LANCZOS_PARAMS_IMMM_H
