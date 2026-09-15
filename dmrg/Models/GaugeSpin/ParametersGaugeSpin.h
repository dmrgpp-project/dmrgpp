// SPDX-FileCopyrightText: Copyright (c) 2009-2012, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 2.0.0]
// [by G.A., Oak Ridge National Laboratory]


/** \ingroup DMRG */
/*@{*/

/*! \file ParametersGaugeSpin.h
 *
 *  Contains the parameters for the Heisenberg model and function
 *  to read them from a file
 *
 */
#ifndef DMRG_PARAMS_MODEL_GAUGESPIN_H
#define DMRG_PARAMS_MODEL_GAUGESPIN_H
#include "ParametersModelBase.h"
#include <PsimagLite/Vector.h>

namespace Dmrg {
//! Heisenberg Model Parameters
template <typename RealType, typename QnType>
struct ParametersGaugeSpin : public ParametersModelBase<RealType, QnType> {

	using BaseType       = ParametersModelBase<RealType, QnType>;
	using VectorRealType = typename PsimagLite::Vector<RealType>::Type;
	// no connectors here, connectors are handled by the geometry
	template <typename IoInputType>
	ParametersGaugeSpin(IoInputType& io)
	    : BaseType(io, false)
	{
		SizeType nsites = 0;
		io.readline(nsites, "TotalNumberOfSites=");

		try {
			magneticFieldV.resize(nsites);
			io.read(magneticFieldV, "MagneticField");
		} catch (std::exception&) { }
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		PsimagLite::String label = label1 + "/ParametersGaugeSpin";
		io.createGroup(label);
		BaseType::write(label, io);
		io.write(label + "/magneticFieldV", magneticFieldV);
	}

	//! Function that prints model parameters to stream os
	friend std::ostream& operator<<(std::ostream& os, const ParametersGaugeSpin& parameters)
	{
		os << "MagneticField=" << parameters.magneticFieldV << "\n";
		os << parameters.targetQuantum;
		return os;
	}

	SizeType       twiceTheSpin;
	VectorRealType magneticFieldV;
};
} // namespace Dmrg

/*@}*/
#endif
