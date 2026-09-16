// SPDX-FileCopyrightText: Copyright (c) 2009-2015, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file ParametersFermionSpinless.h
 *
 *  Contains the parameters for this model
 *
 */
#ifndef DMRG_PARAMS_FERMIONSPINLESS_H
#define DMRG_PARAMS_FERMIONSPINLESS_H
#include "ParametersModelBase.h"

namespace Dmrg {
//! Hubbard Model Parameters
template <typename RealType, typename QnType>
struct ParametersFermionSpinless : public ParametersModelBase<RealType, QnType> {

	using BaseType       = ParametersModelBase<RealType, QnType>;
	using VectorRealType = typename PsimagLite::Vector<RealType>::Type;

	template <typename IoInputType>
	ParametersFermionSpinless(IoInputType& io)
	    : BaseType(io, false)
	{
		SizeType numberOfSites = 0;
		io.readline(numberOfSites, "TotalNumberOfSites=");
		potentialV.resize(numberOfSites);

		io.read(potentialV, "potentialV");
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		PsimagLite::String label = label1 + "/ParametersFermionSpinless";
		io.createGroup(label);
		BaseType::write(label, io);
		io.write(label + "/potentialV", potentialV);
	}

	//! Function that prints model parameters to stream os
	friend std::ostream& operator<<(std::ostream&                    os,
	                                const ParametersFermionSpinless& parameters)
	{
		os << parameters.targetQuantum;
		os << "potentialV\n";
		os << parameters.potentialV;
		return os;
	}

	VectorRealType potentialV;
};
} // namespace Dmrg

/*@}*/
#endif
