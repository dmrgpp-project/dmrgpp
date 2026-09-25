// SPDX-FileCopyrightText: Copyright (c) 2009-2015, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file ParametersHeisenbergAncillaC.h
 *
 *  Contains the parameters for the Heisenberg model and function
 *  to read them from a file
 *
 */
#ifndef DMRG_PARAMS_HEISENBERG_ANCILLAC_H
#define DMRG_PARAMS_HEISENBERG_ANCILLAC_H
#include "ParametersModelBase.h"
#include <PsimagLite/Vector.h>

namespace Dmrg {
//! Heisenberg Model Parameters
template <typename RealType, typename QnType>
struct ParametersHeisenbergAncillaC : public ParametersModelBase<RealType, QnType> {

	using BaseType       = ParametersModelBase<RealType, QnType>;
	using VectorRealType = typename PsimagLite::Vector<RealType>::Type;

	// no connectors here, connectors are handled by the geometry
	template <typename IoInputType>
	ParametersHeisenbergAncillaC(IoInputType& io)
	    : BaseType(io, false)
	{
		io.readline(twiceTheSpin, "HeisenbergTwiceS=");

		try {
			io.read(magneticField, "MagneticField");
		} catch (std::exception&) { }
	}

	template <typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType&, SizeType, PsimagLite::String = "") const
	{
		return 0;
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		PsimagLite::String label = label1 + "/ParametersHeisenbergAncillaC";
		io.createGroup(label);
		BaseType::write(label, io);
		io.write(label + "/twiceTheSpin", twiceTheSpin);
		io.write(label + "/magneticField", magneticField);
	}

	//! Function that prints model parameters to stream os
	friend std::ostream& operator<<(std::ostream&                       os,
	                                const ParametersHeisenbergAncillaC& parameters)
	{
		os << "MagneticField=" << parameters.magneticField << "\n";
		os << "HeisenbergTwiceS=" << parameters.twiceTheSpin << "\n";
		return os;
	}

	SizeType       twiceTheSpin;
	VectorRealType magneticField;
};
} // namespace Dmrg

/*@}*/
#endif
