// SPDX-FileCopyrightText: Copyright (c) 2009-2012-2018, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file ParametersKitaev.h
 *
 *  Contains the parameters for the Kitaev model and function
 *  to read them from a file (started March 2018)
 *
 */
#ifndef PARAMETERS_KITAEV_H
#define PARAMETERS_KITAEV_H
#include "ParametersModelBase.h"
#include <PsimagLite/Vector.h>

namespace Dmrg {
//! Kitaev Model Parameters
// no connectors here, connectors are handled by the geometry
template <typename RealType, typename QnType>
struct ParametersKitaev : public ParametersModelBase<RealType, QnType> {

	using BaseType       = ParametersModelBase<RealType, QnType>;
	using VectorRealType = typename PsimagLite::Vector<RealType>::Type;

	template <typename IoInputType>
	ParametersKitaev(IoInputType& io)
	    : BaseType(io, false)
	{
		SizeType nsites = 0;
		io.readline(nsites, "TotalNumberOfSites=");
		try {
			magneticFieldX.resize(nsites, 0.0);
			io.read(magneticFieldX, "MagneticFieldX");
			std::cerr << "Has MagneticFieldX \n";
		} catch (std::exception&) {
			magneticFieldX.clear();
		}

		try {
			magneticFieldY.resize(nsites, 0.0);
			io.read(magneticFieldY, "MagneticFieldY");
			std::cerr << "Has MagneticFieldY \n";
		} catch (std::exception&) {
			magneticFieldY.clear();
		}

		try {
			magneticFieldZ.resize(nsites, 0.0);
			io.read(magneticFieldZ, "MagneticFieldZ");
			std::cerr << "Has MagneticFieldZ \n";
		} catch (std::exception&) {
			magneticFieldZ.clear();
		}
	}

	template <typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType&, SizeType, PsimagLite::String = "") const
	{
		return 0;
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		PsimagLite::String label = label1 + "/ParametersKitaev";
		io.createGroup(label);
		BaseType::write(label, io);
		io.write(label + "/magneticFieldX", magneticFieldX);
		io.write(label + "/magneticFieldY", magneticFieldY);
		io.write(label + "/magneticFieldZ", magneticFieldZ);
	}

	//! Function that prints model parameters to stream os
	friend std::ostream& operator<<(std::ostream& os, const ParametersKitaev& parameters)
	{
		os << "MagneticFieldX=" << parameters.magneticFieldX << "\n";
		os << "MagneticFieldY=" << parameters.magneticFieldY << "\n";
		os << "MagneticFieldZ=" << parameters.magneticFieldZ << "\n";
		os << parameters.targetQuantum;
		return os;
	}

	VectorRealType magneticFieldX;
	VectorRealType magneticFieldY;
	VectorRealType magneticFieldZ;
};
} // namespace Dmrg

/*@}*/
#endif
