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
#ifndef PARAMS_HUBBARD_HOLSTEIN_SPINLESS_H
#define PARAMS_HUBBARD_HOLSTEIN_SPINLESS_H
#include "ParametersModelBase.h"
#include <PsimagLite/Matrix.h>
#include <stdexcept>
#include <vector>

namespace Dmrg {
//! FeAs Model Parameters
template <typename ComplexOrRealType, typename QnType>
struct ParametersHubbardHolsteinSpinless : public ParametersModelBase<ComplexOrRealType, QnType> {
	// no connections here please!!
	// connections are handled by the geometry

	using RealType       = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using BaseType       = ParametersModelBase<ComplexOrRealType, QnType>;
	using VectorRealType = typename PsimagLite::Vector<RealType>::Type;
	using VectorType = typename PsimagLite::Vector<PsimagLite::Matrix<ComplexOrRealType>>::Type;

	template <typename IoInputType>
	ParametersHubbardHolsteinSpinless(IoInputType& io)
	    : BaseType(io, false)
	    , oStruncPhonons(0)
	    , oStruncSite(0)
	{

		SizeType nsites = 0;
		io.readline(nsites, "TotalNumberOfSites=");
		PsimagLite::String model;
		io.readline(model, "Model=");
		potentialFV.resize(nsites, 0.0);
		potentialPV.resize(nsites, 0.0);
		io.readline(numberphonons, "NumberPhonons=");
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

		potentialDV.resize(nsites, 0.0);
		try {
			io.read(potentialDV, "potentialDV");
		} catch (...) {
			potentialDV.clear();
		}

		try {
			io.readline(oStruncPhonons, "OneSiteTruncationPhononsMax=");
			io.readline(oStruncSite, "OneSiteTruncationSite=");
		} catch (...) { }

		if (oStruncPhonons > 0 && oStruncSite == 0)
			err("OneSiteTruncationSite cannot be zero\n");

		if (model == "HubbardHolsteinSpinless" || model == "HubbardHolsteinSpinlessSSH") {
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
		PsimagLite::String label = label1 + "/ParametersHubbardHolsteinSpinless";
		io.createGroup(label);
		BaseType::write(label, io);
		io.write(label + "/numberphonons", numberphonons);
		io.write(label + "/lambdaFP", lambdaFP);
		io.write(label + "/potentialFV", potentialFV);
		io.write(label + "/potentialPV", potentialPV);
		io.write(label + "/potentialDV", potentialDV);
	}

	//! Function that prints model parameters to stream os
	friend std::ostream& operator<<(std::ostream&                            os,
	                                const ParametersHubbardHolsteinSpinless& parameters)
	{
		os << "NumberPhonons=" << parameters.numberphonons << "\n";
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
	VectorRealType lambdaFP;
	// Onsite potential values
	VectorRealType potentialFV;
	VectorRealType potentialPV;
	VectorRealType potentialDV;
};
} // namespace Dmrg

/*@}*/
#endif
