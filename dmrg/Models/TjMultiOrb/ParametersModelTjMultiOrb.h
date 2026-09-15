// SPDX-FileCopyrightText: Copyright (c) 2009-2012, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 2.0.0]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file ParametersModelTjMultiOrb.h
 *
 *  Contains the parameters for the Hubbard model and function to
 *  read them from a file
 *
 */
#ifndef DMRG_PARAMS_TJMULTIORB_H
#define DMRG_PARAMS_TJMULTIORB_H
#include "ParametersModelBase.h"

namespace Dmrg {
//! ModelTjMultiOrb Parameters
template <typename RealType, typename QnType>
struct ParametersModelTjMultiOrb : public ParametersModelBase<RealType, QnType> {

	using BaseType = ParametersModelBase<RealType, QnType>;

	template <typename IoInputType>
	ParametersModelTjMultiOrb(IoInputType& io)
	    : BaseType(io, false)
	    , reinterpretAndTruncate(0)
	{
		io.read(potentialV, "potentialV");
		io.readline(orbitals, "Orbitals=");

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
		PsimagLite::String label = label1 + "/ParametersModelTjMultiOrb";
		io.createGroup(label);
		BaseType::write(label, io);
		io.write(label + "/potentialV", potentialV);
		io.write(label + "/orbitals", orbitals);
		io.write(label + "/reinterpretAndTruncate", reinterpretAndTruncate);
	}

	//! Function that prints model parameters to stream os
	friend std::ostream& operator<<(std::ostream&                    os,
	                                const ParametersModelTjMultiOrb& parameters)
	{
		os << "potentialV\n";
		os << parameters.potentialV;
		os << "orbitals=" << parameters.orbitals << "\n";
		os << "JHundInfinity=" << parameters.reinterpretAndTruncate << "\n";
		return os;
	}

	// Do not include here connection parameters
	typename PsimagLite::Vector<RealType>::Type potentialV;
	SizeType                                    orbitals;
	SizeType                                    reinterpretAndTruncate;
};
} // namespace Dmrg

/*@}*/
#endif
