// SPDX-FileCopyrightText: Copyright (c) 2009-2018, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file ParametersHubbardMultiBand.h
 *
 *  Contains the parameters for the HubbardMultiBand model
 *
 */
#ifndef PARAMS_HUBBARD_MULTI_BAND_H
#define PARAMS_HUBBARD_MULTI_BAND_H
#include "ParametersModelBase.h"
#include <PsimagLite/Matrix.h>
#include <stdexcept>
#include <vector>

namespace Dmrg {
//! FeAs Model Parameters
template <typename ComplexOrRealType, typename QnType>
struct ParametersHubbardMultiBand : public ParametersModelBase<ComplexOrRealType, QnType> {
	// no connections here please!!
	// connections are handled by the geometry

	using BaseType       = ParametersModelBase<ComplexOrRealType, QnType>;
	using RealType       = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using VectorRealType = typename PsimagLite::Vector<RealType>::Type;
	using VectorType = typename PsimagLite::Vector<PsimagLite::Matrix<ComplexOrRealType>>::Type;

	template <typename IoInputType>
	ParametersHubbardMultiBand(IoInputType& io)
	    : BaseType(io, false)
	{
		io.readline(orbitals, "Orbitals=");
		io.read(hubbardU, "hubbardU");
		io.read(potentialV, "potentialV");
		SizeType h = 1;
		try {
			io.readline(h, "NumberOfHoppingOrbitalMatrices=");
		} catch (std::exception&) { }

		SizeType sites = 0;
		io.readline(sites, "TotalNumberOfSites=");
		if (h != 1 && h != sites)
			err("NumberOfHoppingOrbitalMatrices=1 or =numberOfSites\n");

		for (SizeType i = 0; i < h; ++i) {
			PsimagLite::Matrix<ComplexOrRealType> m;
			io.read(m, "hopOnSite");
			hopOnSite.push_back(m);
		}
	}

	template <typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType&, SizeType, PsimagLite::String = "") const
	{
		return 0;
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		PsimagLite::String label = label1 + "/ParametersHubbardMultiBand";
		io.createGroup(label);
		BaseType::write(label, io);
		io.write(label + "/orbitals", orbitals);
		io.write(label + "/hubbardU", hubbardU);
		io.write(label + "/potentialV", potentialV);
		io.write(label + "/hopOnSite", hopOnSite);
	}

	//! Function that prints model parameters to stream os
	friend std::ostream& operator<<(std::ostream&                     os,
	                                const ParametersHubbardMultiBand& parameters)
	{
		os << "Orbitals=" << parameters.orbitals << "\n";
		os << "hubbardU\n";
		os << parameters.hubbardU;
		os << "potentialV\n";
		os << parameters.potentialV;

		return os;
	}

	SizeType       orbitals;
	VectorRealType hubbardU;
	// Onsite potential values, one for each site
	VectorRealType potentialV;
	VectorType     hopOnSite;
};
} // namespace Dmrg

/*@}*/
#endif
