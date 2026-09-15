// SPDX-FileCopyrightText: Copyright (c) 2009-2014, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file MettsParams.h
 *
 *  This is a structure to represent the parameters of the TimeStep Evolution
 *  algorithm. Don't add functions to this class because
 *  this class's data is all public
 */
#ifndef METTS_PARAMS_H
#define METTS_PARAMS_H

#include "TargetParamsTimeVectors.h"

namespace Dmrg {

// Coordinates reading of TargetSTructure from input file
template <typename ModelType> class MettsParams : public TargetParamsTimeVectors<ModelType> {

	using BaseType       = TargetParamsTimeVectors<ModelType>;
	using VectorSizeType = typename BaseType::VectorSizeType;

public:

	using RealType             = typename ModelType::RealType;
	using OperatorType         = typename ModelType::OperatorType;
	using TimeVectorParamsType = TargetParamsTimeVectors<ModelType>;
	using SparseMatrixType     = typename OperatorType::StorageType;

	template <typename IoInputter>
	MettsParams(IoInputter& io, PsimagLite::String targeting, const ModelType& model)
	    : TimeVectorParamsType(io, targeting, model)
	{
		io.readline(beta, "BetaDividedByTwo=");
		io.readline(rngSeed, "TSPRngSeed=");
		io.readline(collapse, "MettsCollapse=");
		try {
			io.read(pure, "MettsPure");
		} catch (std::exception& e) { }

		SizeType n = model.superGeometry().numberOfSites();
		if (pure.size() > 0 && pure.size() != n) {
			PsimagLite::String msg("MettsParams: If provided, MettsPure must be");
			msg += " a vector of " + ttos(n) + " entries.\n";
			throw PsimagLite::RuntimeError(msg);
		}

		this->noOperator(false);
	}

	int long           rngSeed;
	RealType           beta;
	PsimagLite::String collapse;
	VectorSizeType     pure;
}; // class MettsParams

template <typename ModelType>
inline std::ostream& operator<<(std::ostream& os, const MettsParams<ModelType>& t)
{
	os << "TargetParams.type=Metts";
	const typename MettsParams<ModelType>::TimeVectorParamsType& tp = t;
	os << tp;
	os << "BetaDividedByTwo=" << t.beta << "\n";
	os << "TSPRngSeed=" << t.rngSeed << "\n";
	os << "MettsCollapse=" << t.collapse << "\n";
	os << "MettsPure=" << t.pure << "\n";
	return os;
}

} // namespace Dmrg

/*@}*/
#endif // METTS_PARAMS_H
