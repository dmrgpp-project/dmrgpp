// SPDX-FileCopyrightText: Copyright (c) 2009-2014, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 2.0.0]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file ParametersModelHubbard.h
 *
 *  Contains the parameters for the Hubbard model and function to read them from a file
 *
 */
#ifndef PARAMETERSMODELHUBBARD_H
#define PARAMETERSMODELHUBBARD_H
#include "InputCheck.h"
#include "ParametersModelBase.h"
#include <PsimagLite/InputNg.h>

namespace Dmrg {
//! Hubbard Model Parameters
template <typename RealType, typename QnType>
struct ParametersModelHubbard : public ParametersModelBase<RealType, QnType> {

	using BaseType         = ParametersModelBase<RealType, QnType>;
	using IoInputType      = PsimagLite::InputNg<InputCheck>::Readable;
	using VectorStringType = PsimagLite::Vector<PsimagLite::String>::Type;
	using VectorRealType   = std::vector<RealType>;

	ParametersModelHubbard(IoInputType& io)
	    : BaseType(io, false)
	{
		SizeType nsites = 0;
		io.readline(nsites, "TotalNumberOfSites=");
		hubbardU.resize(nsites, 0.0);
		potentialV.resize(2 * nsites, 0.0);
		io.read(hubbardU, "hubbardU");
		io.read(potentialV, "potentialV");
		try {
			anisotropy.resize(nsites, 0.0);
			io.read(anisotropy, "AnisotropyD");
			std::cerr << "Has AnisotropyD\n";
		} catch (std::exception&) {
			anisotropy.clear();
		}

		if (anisotropy.size() > 0 && anisotropy.size() != nsites) {
			throw PsimagLite::RuntimeError("AnisotropyD size must be " + ttos(nsites)
			                               + " entries long, if provided.\n");
		}

		try {
			magneticX.resize(nsites, 0.0);
			io.read(magneticX, "MagneticFieldX");
			std::cerr << "Has MagneticFieldX\n";
		} catch (std::exception&) {
			magneticX.clear();
		}

		if (magneticX.size() > 0 && magneticX.size() != nsites) {
			throw PsimagLite::RuntimeError("MagneticFieldX size must be " + ttos(nsites)
			                               + " entries long, if provided.\n");
		}

		onSiteHaddLegacy = readOldT(io, nsites);
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		PsimagLite::String label = label1 + "/ParametersModelHubbard";
		io.createGroup(label);
		BaseType::write(label, io);
		io.write(label + "/hubbardU", hubbardU);
		io.write(label + "/potentialV", potentialV);
		io.write(label + "/anisotropy", anisotropy);
		io.write(label + "/magneticX", magneticX);
	}

	VectorStringType readOldT(IoInputType& io, SizeType nsites)
	{
		VectorRealType potentialTlegacy;
		try {
			io.read(potentialTlegacy, "PotentialT");
			std::cerr << "Has PotentialT\n";
		} catch (std::exception&) {
			return VectorStringType();
		}

		RealType omega = 0;
		try {
			io.readline(omega, "omega=");
		} catch (std::exception&) { }

		RealType phase = 0;
		try {
			io.readline(phase, "phase=");
		} catch (std::exception&) { }

		// c means cosine below
		const PsimagLite::String function
		    = "*(c:+:*:%t:" + ttos(omega) + ":" + ttos(phase) + ")*";
		const PsimagLite::String nup   = function + "nup";
		const PsimagLite::String ndown = function + "ndown";

		VectorStringType potentialTv(nsites);
		for (SizeType site = 0; site < nsites; ++site) {
			const RealType           val = potentialTlegacy[site];
			const PsimagLite::String plusSignOrNot
			    = ((val < 0) && (site > 0)) ? "+" : "";
			PsimagLite::String expression = ttos(val) + nup + " + ";
			expression += plusSignOrNot + ttos(val) + ndown;
			potentialTv[site] = ProgramGlobals::killSpaces(expression);
		}

		return potentialTv;
	}

	VectorRealType hubbardU;
	VectorRealType potentialV;
	VectorRealType anisotropy;
	VectorRealType magneticX;

	// for time-dependent H:
	VectorStringType onSiteHaddLegacy;
};
} // namespace Dmrg

/*@}*/
#endif
