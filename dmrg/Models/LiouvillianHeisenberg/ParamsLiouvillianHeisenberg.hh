// SPDX-FileCopyrightText: Copyright (c) 2009, 2012-2026, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 6+]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file ParamsLiouvillianHeisenberg.hh
 *
 *  Contains the parameters for the LiouvillianHeisenberg model and functions
 *  to read them from a file
 *
 */
#ifndef DMRG_PARAMS_LIOUVILLIANHEISENBERG_H
#define DMRG_PARAMS_LIOUVILLIANHEISENBERG_H
#include "ParametersModelBase.h"
#include "ProgramGlobals.h"
#include <PsimagLite/Vector.h>

namespace Dmrg {
//! Heisenberg Model Parameters
template <typename RealType, typename QnType>
struct ParamsLiouvillianHeisenberg : public ParametersModelBase<RealType, QnType> {

	using BaseType       = ParametersModelBase<RealType, QnType>;
	using VectorRealType = typename PsimagLite::Vector<RealType>::Type;

	enum class Ancillas
	{
		TOGETHER,
		SEPARATE
	};

	// no connectors here, connectors are handled by the geometry
	template <typename IoInputType>
	ParamsLiouvillianHeisenberg(IoInputType& io)
	    : BaseType(io, false)
	{
		PsimagLite::String model;
		io.readline(model, "Model=");

		io.readline(twiceTheSpin, "HeisenbergTwiceS=");

		{
			std::string tmp;
			io.readline(tmp, "Ancillas");
			tmp = ProgramGlobals::toLower(tmp);
			if (tmp == "together") {
				ancillas = Ancillas::TOGETHER;
			} else if (tmp == "separate") {
				ancillas = Ancillas::SEPARATE;
			} else {
				err("Ancillas= together or separate but not " + tmp + "\n");
			}
		}

		SizeType nsites = 0;
		io.readline(nsites, "TotalNumberOfSites=");

		bath_gamma.resize(nsites);
		io.read(bath_gamma, "BathGamma");

		bath_f.resize(nsites);
		io.read(bath_f, "BathF");

		try {
			magneticFieldZ.resize(nsites);
			io.read(magneticFieldZ, "MagneticFieldZ");
		} catch (std::exception&) {
			magneticFieldZ.clear();
		}
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		PsimagLite::String label = label1 + "/ParamsLiouvillianHeisenberg";
		io.createGroup(label);
		BaseType::write(label, io);
		io.write(label + "/twiceTheSpin", twiceTheSpin);
		io.write(label + "/magneticFieldZ", magneticFieldZ);
	}

	static void checkMagneticField(SizeType s, unsigned char c, SizeType n)
	{
		if (s == 0 || s == n)
			return;

		PsimagLite::String msg("LiouvillianHeisenberg: If provided, ");
		msg += " MagneticField" + ttos(c) + " must be a vector of " + ttos(n)
		    + " entries.\n";
		err(msg);
	}

	//! Function that prints model parameters to stream os
	friend std::ostream& operator<<(std::ostream&                      os,
	                                const ParamsLiouvillianHeisenberg& parameters)
	{
		if (!parameters.magneticFieldZ.empty())
			os << "MagneticFieldZ=" << parameters.magneticFieldZ << "\n";

		os << "HeisenbergTwiceS=" << parameters.twiceTheSpin << "\n";
		os << parameters.targetQuantum;
		return os;
	}

	SizeType       twiceTheSpin;
	VectorRealType magneticFieldZ;
	VectorRealType bath_gamma;
	VectorRealType bath_f;
	Ancillas       ancillas = Ancillas::TOGETHER;
};
} // namespace Dmrg

/*@}*/
#endif
