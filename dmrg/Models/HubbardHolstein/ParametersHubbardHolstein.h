// SPDX-FileCopyrightText: Copyright (c) 2009-2015, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file ParametersHubbardHolstein.h
 *
 *
 */
#ifndef PARAMS_HUBBARD_HOLSTEIN_H
#define PARAMS_HUBBARD_HOLSTEIN_H
#include "ParametersModelBase.h"
#include <PsimagLite/Matrix.h>
#include <stdexcept>
#include <vector>

namespace Dmrg {
//! FeAs Model Parameters
template <typename ComplexOrRealType, typename QnType>
struct ParametersHubbardHolstein : public ParametersModelBase<ComplexOrRealType, QnType> {
	// no connections here please!!
	// connections are handled by the geometry

	using RealType       = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using BaseType       = ParametersModelBase<ComplexOrRealType, QnType>;
	using VectorRealType = typename PsimagLite::Vector<RealType>::Type;
	using VectorType = typename PsimagLite::Vector<PsimagLite::Matrix<ComplexOrRealType>>::Type;

	template <typename IoInputType>
	ParametersHubbardHolstein(IoInputType& io)
	    : BaseType(io, false)
	    , oStruncPhonons(0)
	    , oStruncSite(0)
	{

		SizeType nsites = 0;
		io.readline(nsites, "TotalNumberOfSites=");
		PsimagLite::String model;
		io.readline(model, "Model=");
		hubbardFU.resize(nsites, 0.0);
		potentialFV.resize(2 * nsites, 0.0);
		potentialPV.resize(nsites, 0.0);
		io.readline(numberphonons, "NumberPhonons=");
		io.read(hubbardFU, "hubbardFU");
		io.read(potentialFV, "potentialFV");
		io.read(potentialPV, "potentialPV");

		lambdaFP.resize(nsites, 0.0);
		bool hasLambdaFP = false;
		try {
			io.read(lambdaFP, "lambdaFP");
			hasLambdaFP = true;
		} catch (...) {
			lambdaFP.clear();
		}

		try {
			io.readline(oStruncPhonons, "OneSiteTruncationPhononsMax=");
			io.readline(oStruncSite, "OneSiteTruncationSite=");
		} catch (...) { }

		if (oStruncPhonons > 0 && oStruncSite == 0)
			err("OneSiteTruncationSite cannot be zero\n");

		if (model == "HubbardHolstein" || model == "HubbardHolsteinSSH") {
			if (!hasLambdaFP)
				err("HubbardHolstein: must have lambdaFP vector in input file\n");
		} else {
			if (hasLambdaFP)
				err("HolsteinThin: lambdaFP should be given as a connection\n");
		}
	}

	template <typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType&, SizeType, PsimagLite::String = "") const
	{
		return 0;
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		PsimagLite::String label = label1 + "/ParametersHubbardHolstein";
		io.createGroup(label);
		BaseType::write(label, io);
		io.write(label + "/numberphonons", numberphonons);
		io.write(label + "/hubbardFU", hubbardFU);
		io.write(label + "/lambdaFP", lambdaFP);
		io.write(label + "/potentialFV", potentialFV);
		io.write(label + "/potentialPV", potentialPV);
	}

	//! Function that prints model parameters to stream os
	friend std::ostream& operator<<(std::ostream&                    os,
	                                const ParametersHubbardHolstein& parameters)
	{
		os << "NumberPhonons=" << parameters.numberphonons << "\n";

		os << "hubbardFU\n";
		os << parameters.hubbardFU;
		os << "lambdaFP\n";
		os << parameters.lambdaFP;

		os << "potentialFV\n";
		os << parameters.potentialFV;
		os << "potentialPV\n";
		os << parameters.potentialPV;
		return os;
	}

	SizeType       numberphonons;
	SizeType       oStruncPhonons;
	SizeType       oStruncSite;
	VectorRealType hubbardFU;
	VectorRealType lambdaFP;
	// Onsite potential values
	VectorRealType potentialFV;
	VectorRealType potentialPV;
};
} // namespace Dmrg

/*@}*/
#endif
