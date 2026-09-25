// SPDX-FileCopyrightText: Copyright (c) 2009,-2012 UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 2.0.0]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/
/** \file Parallel2PointCorrelations.h
 */

#ifndef PARALLEL_2POINT_CORRELATIONS_H
#define PARALLEL_2POINT_CORRELATIONS_H

#include "ProgramGlobals.h"
#include <PsimagLite/Concurrency.h>
#include <PsimagLite/GetBraOrKet.h>
#include <PsimagLite/Matrix.h>
#include <PsimagLite/NotMpi.h>

namespace Dmrg {

template <typename TwoPointCorrelationsType> class Parallel2PointCorrelations {

public:

	using MatrixType       = typename TwoPointCorrelationsType::MatrixType;
	using SparseMatrixType = typename TwoPointCorrelationsType::SparseMatrixType;
	using FieldType        = typename MatrixType::value_type;
	using ConcurrencyType  = PsimagLite::Concurrency;
	using PairType         = std::pair<SizeType, SizeType>;
	using RealType         = typename PsimagLite::Real<FieldType>::Type;
	using BraketType       = typename TwoPointCorrelationsType::BraketType;

	Parallel2PointCorrelations(MatrixType&                                        w,
	                           const TwoPointCorrelationsType&                    twopoint,
	                           const typename PsimagLite::Vector<PairType>::Type& pairs,
	                           const BraketType&                                  braket,
	                           ProgramGlobals::FermionOrBosonEnum                 fermionicSign,
	                           const PsimagLite::GetBraOrKet&                     bra,
	                           const PsimagLite::GetBraOrKet&                     ket)
	    : w_(w)
	    , twopoint_(twopoint)
	    , pairs_(pairs)
	    , braket_(braket)
	    , fermionicSign_(fermionicSign)
	    , bra_(bra)
	    , ket_(ket)
	{ }

	void doTask(SizeType taskNumber, SizeType)
	{
		SizeType i = pairs_[taskNumber].first;
		SizeType j = pairs_[taskNumber].second;
		w_(i, j)   = twopoint_.calcCorrelation(i, j, braket_, fermionicSign_, bra_, ket_);
	}

	SizeType tasks() const { return pairs_.size(); }

private:

	MatrixType&                                        w_;
	const TwoPointCorrelationsType&                    twopoint_;
	const typename PsimagLite::Vector<PairType>::Type& pairs_;
	const BraketType&                                  braket_;
	const ProgramGlobals::FermionOrBosonEnum           fermionicSign_;
	const PsimagLite::GetBraOrKet&                     bra_;
	const PsimagLite::GetBraOrKet&                     ket_;
}; // class Parallel2PointCorrelations
} // namespace Dmrg

/*@}*/
#endif // PARALLEL_2POINT_CORRELATIONS_H
