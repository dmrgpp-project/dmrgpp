// SPDX-FileCopyrightText: Copyright (c) 2009-2014, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file ParametersImmm.h
 *
 *  Contains the parameters for the Immm model and function to read them
 *
 */
#ifndef PARAMETERS_IMMM_H
#define PARAMETERS_IMMM_H
#include "ParametersModelBase.h"

namespace Dmrg {
template <typename RealType, typename QnType>
struct ParametersImmm : public ParametersModelBase<RealType, QnType> {
	// no connections here please!!
	// connections are handled by the geometry

	using BaseType = ParametersModelBase<RealType, QnType>;

	template <typename IoInputType>
	ParametersImmm(IoInputType& io)
	    : BaseType(io, false)
	{

		io.read(hubbardU, "hubbardU");
		io.read(potentialV, "potentialV");
		io.readline(minOxygenElectrons, "MinOxygenElectrons=");
	}

	template <typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType& mres, SizeType, PsimagLite::String msg = "") const
	{
		PsimagLite::String str = msg;
		str += "ParametersImmm";

		const char* start = reinterpret_cast<const char*>(this);
		const char* end   = reinterpret_cast<const char*>(&potentialV);
		SizeType    total = mres.memResolv(&hubbardU, end - start, str + " hubbardU");

		start = end;
		end   = reinterpret_cast<const char*>(&minOxygenElectrons);
		total += mres.memResolv(&potentialV, end - start, str + " potentialV");

		total += mres.memResolv(
		    &minOxygenElectrons, sizeof(*this) - total, str + " minOxygenElectrons");

		return total;
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		PsimagLite::String label = label1 + "/ParametersImmm";
		io.createGroup(label);
		BaseType::write(label, io);
		io.write(label + "/hubbardU", hubbardU);
		io.write(label + "/potentialV", potentialV);
		io.write(label + "/minOxygenElectrons", minOxygenElectrons);
	}

	//! Function that prints model parameters to stream os
	friend std::ostream& operator<<(std::ostream& os, const ParametersImmm& parameters)
	{
		os << "hubbardU\n";
		os << parameters.hubbardU;
		os << "potentialV\n";
		os << parameters.potentialV;
		os << "MinOxygenElectrons=" << parameters.minOxygenElectrons << "\n";
		return os;
	}

	typename PsimagLite::Vector<RealType>::Type hubbardU;
	typename PsimagLite::Vector<RealType>::Type potentialV;
	SizeType                                    minOxygenElectrons;
};
} // namespace Dmrg

/*@}*/
#endif // PARAMETERS_IMMM_H
