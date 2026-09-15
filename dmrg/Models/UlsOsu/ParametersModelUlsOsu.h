// SPDX-FileCopyrightText: Copyright (c) 2009-2012, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 2.0.0]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file ParametersModelUlsOsu.h
 *
 *  Contains the parameters for the Hubbard model and function to
 *  read them from a file
 *
 */
#ifndef DMRG_PARAMS_UlsOsu_H
#define DMRG_PARAMS_UlsOsu_H
#include "ParametersModelBase.h"

namespace Dmrg {
//! ModelUlsOsu Parameters
template <typename RealType, typename QnType>
struct ParametersModelUlsOsu : public ParametersModelBase<RealType, QnType> {

	using BaseType       = ParametersModelBase<RealType, QnType>;
	using VectorRealType = typename PsimagLite::Vector<RealType>::Type;

	template <typename IoInputType>
	ParametersModelUlsOsu(IoInputType& io)
	    : BaseType(io, false)
	{
		SizeType nsites = 0;
		io.readline(nsites, "TotalNumberOfSites=");
		io.readline(orbitals, "Orbitals=");

		magneticFieldX.resize(nsites, 0.0);
		magneticFieldY.resize(nsites, 0.0);
		magneticFieldZ.resize(nsites, 0.0);

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
	}

	template <typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType&, SizeType, PsimagLite::String = "") const
	{
		return 0;
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		PsimagLite::String label = label1 + "/ParametersModelUlsOsu";
		io.createGroup(label);
		BaseType::write(label, io);
		io.write(label + "/magneticFieldX", magneticFieldX);
		io.write(label + "/magneticFieldY", magneticFieldY);
		io.write(label + "/magneticFieldZ", magneticFieldZ);
	}

	//! Function that prints model parameters to stream os
	friend std::ostream& operator<<(std::ostream& os, const ParametersModelUlsOsu& parameters)
	{
		os << "MagneticFieldX=" << parameters.magneticFieldX << "\n";
		os << "MagneticFieldY=" << parameters.magneticFieldY << "\n";
		os << "MagneticFieldZ=" << parameters.magneticFieldZ << "\n";
		return os;
	}

	// Do not include here connection parameters
	SizeType orbitals;

	VectorRealType magneticFieldX;
	VectorRealType magneticFieldY;
	VectorRealType magneticFieldZ;
};
} // namespace Dmrg

/*@}*/
#endif
