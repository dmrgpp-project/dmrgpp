// SPDX-FileCopyrightText: Copyright (c) 2009-2012, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 2.0.0]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file ParametersModelTjAnisotropic.h
 *
 *  Contains the parameters for the Hubbard model and function to
 *  read them from a file
 *
 */
#ifndef DMRG_PARAMS_TJAnisotropic_H
#define DMRG_PARAMS_TJAnisotropic_H
#include "ParametersModelBase.h"

namespace Dmrg {
//! ModelTjAnisotropic Parameters
template <typename RealType, typename QnType>
struct ParametersModelTjAnisotropic : public ParametersModelBase<RealType, QnType> {

	using BaseType       = ParametersModelBase<RealType, QnType>;
	using VectorRealType = typename PsimagLite::Vector<RealType>::Type;

	template <typename IoInputType>
	ParametersModelTjAnisotropic(IoInputType& io)
	    : BaseType(io, false)
	    , reinterpretAndTruncate(0)
	{
		SizeType nsites = 0;
		io.readline(nsites, "TotalNumberOfSites=");
		io.readline(orbitals, "Orbitals=");

		potentialV.resize(nsites * 2, 0.0);
		magneticFieldX.resize(nsites, 0.0);
		magneticFieldY.resize(nsites, 0.0);
		magneticFieldZ.resize(nsites, 0.0);

		try {
			io.read(potentialV, "potentialV");
			std::cerr << "Has potentialV \n";
		} catch (std::exception&) { }

		try {
			io.read(magneticFieldX, "MagneticFieldX");
			std::cerr << "Has MagneticFieldX \n";
		} catch (std::exception&) { }

		try {
			io.read(magneticFieldY, "MagneticFieldY");
			std::cerr << "Has MagneticFieldY \n";
		} catch (std::exception&) { }

		try {
			io.read(magneticFieldZ, "MagneticFieldZ");
			std::cerr << "Has MagneticFieldZ \n";
		} catch (std::exception&) { }

		try {
			io.readline(reinterpretAndTruncate, "JHundInfinity=");
		} catch (std::exception&) { }

		if (orbitals != 2 && reinterpretAndTruncate > 0)
			throw PsimagLite::RuntimeError(
			    "JHundInfinity>0 only possible for orbitals==2\n");

		if (reinterpretAndTruncate > 3)
			throw PsimagLite::RuntimeError(
			    "JHundInfinity must be less or equal to 3\n");
	}

	template <typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType&, SizeType, PsimagLite::String = "") const
	{
		return 0;
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		PsimagLite::String label = label1 + "/ParametersModelTjAnisotropic";
		io.createGroup(label);
		BaseType::write(label, io);
		io.write(label + "/potentialV", potentialV);
		io.write(label + "/orbitals", orbitals);
		io.write(label + "/reinterpretAndTruncate", reinterpretAndTruncate);
		io.write(label + "/magneticFieldX", magneticFieldX);
		io.write(label + "/magneticFieldY", magneticFieldY);
		io.write(label + "/magneticFieldZ", magneticFieldZ);
	}

	//! Function that prints model parameters to stream os
	friend std::ostream& operator<<(std::ostream&                       os,
	                                const ParametersModelTjAnisotropic& parameters)
	{
		os << "potentialV\n";
		os << parameters.potentialV;
		os << "orbitals=" << parameters.orbitals << "\n";
		os << "JHundInfinity=" << parameters.reinterpretAndTruncate << "\n";
		os << "MagneticFieldX=" << parameters.magneticFieldX << "\n";
		os << "MagneticFieldY=" << parameters.magneticFieldY << "\n";
		os << "MagneticFieldZ=" << parameters.magneticFieldZ << "\n";
		return os;
	}

	// Do not include here connection parameters
	VectorRealType potentialV;
	VectorRealType magneticFieldX;
	VectorRealType magneticFieldY;
	VectorRealType magneticFieldZ;

	SizeType orbitals;
	SizeType reinterpretAndTruncate;
};
} // namespace Dmrg

/*@}*/
#endif
