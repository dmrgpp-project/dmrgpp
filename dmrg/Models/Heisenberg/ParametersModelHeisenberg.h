// SPDX-FileCopyrightText: Copyright (c) 2009-2012, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 2.0.0]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file ParametersModelHeisenberg.h
 *
 *  Contains the parameters for the Heisenberg model and function
 *  to read them from a file
 *
 */
#ifndef PARAMETERSMODELHEISENBERG_H
#define PARAMETERSMODELHEISENBERG_H
#include "../../Engine/ParametersModelBase.h"
#include <PsimagLite/Vector.h>

namespace Dmrg {
//! Heisenberg Model Parameters
template <typename RealType, typename QnType>
struct ParametersModelHeisenberg : public ParametersModelBase<RealType, QnType> {

	using BaseType       = ParametersModelBase<RealType, QnType>;
	using VectorRealType = typename PsimagLite::Vector<RealType>::Type;
	// no connectors here, connectors are handled by the geometry
	template <typename IoInputType>
	ParametersModelHeisenberg(IoInputType& io)
	    : BaseType(io, false)
	    , twiceTheSpinBorder(0)
	{
		PsimagLite::String model;
		io.readline(model, "Model=");

		io.readline(twiceTheSpin, "HeisenbergTwiceS=");

		if (model == "HeisenbergMix")
			io.readline(twiceTheSpinBorder, "HeisenbergTwiceSborder=");

		SizeType nsites = 0;
		io.readline(nsites, "TotalNumberOfSites=");

		try {
			magneticFieldX.resize(nsites);
			io.read(magneticFieldX, "MagneticFieldX");
		} catch (std::exception&) {
			magneticFieldX.clear();
		}

		try {
			magneticFieldZ.resize(nsites);
			io.read(magneticFieldZ, "MagneticFieldZ");
		} catch (std::exception&) {
			magneticFieldZ.clear();
		}

		// throw if supplying MagneticField label
		bool invalidLabel = false;
		try {
			VectorRealType tmpVector;
			io.read(tmpVector, "MagneticField=");
			invalidLabel = true;
		} catch (std::exception&) { }

		if (invalidLabel) {
			throw PsimagLite::RuntimeError(
			    "MagneticField label is no longer supported.\n"
			    + PsimagLite::String("Please use MagneticField[XZ] instead\n"));
		}

		// throw if supplying MagneticFieldDirection label
		try {
			PsimagLite::String tmpStr;
			io.readline(tmpStr, "MagneticFieldDirection=");
			invalidLabel = true;
		} catch (std::exception&) { }

		if (invalidLabel) {
			throw PsimagLite::RuntimeError(
			    "MagneticFieldDirection label is no longer supported.\n"
			    + PsimagLite::String("Please use MagneticField[XZ] instead\n"));
		}

		try {
			io.read(anisotropyD, "AnisotropyD");
		} catch (std::exception&) { }

		try {
			io.read(anisotropyE, "AnisotropyE");
		} catch (std::exception&) { }
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		PsimagLite::String label = label1 + "/ParametersModelHeisenberg";
		io.createGroup(label);
		BaseType::write(label, io);
		io.write(label + "/twiceTheSpin", twiceTheSpin);
		io.write(label + "/magneticFieldX", magneticFieldX);
		io.write(label + "/magneticFieldZ", magneticFieldZ);
		io.write(label + "/anisotropyD", anisotropyD);
		io.write(label + "/anisotropyE", anisotropyE);
	}

	static void checkMagneticField(SizeType s, unsigned char c, SizeType n)
	{
		if (s == 0 || s == n)
			return;

		PsimagLite::String msg("ModelHeisenberg: If provided, ");
		msg += " MagneticField" + ttos(c) + " must be a vector of " + ttos(n)
		    + " entries.\n";
		err(msg);
	}

	//! Function that prints model parameters to stream os
	friend std::ostream& operator<<(std::ostream&                    os,
	                                const ParametersModelHeisenberg& parameters)
	{
		if (!parameters.magneticFieldX.empty())
			os << "MagneticFieldX=" << parameters.magneticFieldX << "\n";
		if (!parameters.magneticFieldZ.empty())
			os << "MagneticFieldZ=" << parameters.magneticFieldZ << "\n";

		os << "AnisotropyD=" << parameters.anisotropy << "\n";
		os << "HeisenbergTwiceS=" << parameters.twiceTheSpin << "\n";
		os << parameters.targetQuantum;
		return os;
	}

	SizeType       twiceTheSpin;
	SizeType       twiceTheSpinBorder;
	VectorRealType magneticFieldX;
	VectorRealType magneticFieldZ;
	VectorRealType anisotropyD;
	VectorRealType anisotropyE;
};
} // namespace Dmrg

/*@}*/
#endif
