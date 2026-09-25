// SPDX-FileCopyrightText: Copyright (c) 2013, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5..0]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file MultiPointCorrelations.h
 *
 *
 *
 */
#ifndef MULTI_POINT_CORRELATIONS_H
#define MULTI_POINT_CORRELATIONS_H
#include "VectorWithOffsets.h" // for operator*
#include <PsimagLite/CrsMatrix.h>

namespace Dmrg {

template <typename CorrelationsSkeletonType> class MultiPointCorrelations {

	using ObserverHelperType     = typename CorrelationsSkeletonType::ObserverHelperType;
	using VectorType             = typename ObserverHelperType::VectorType;
	using VectorWithOffsetType   = typename ObserverHelperType::VectorWithOffsetType;
	using BasisWithOperatorsType = typename ObserverHelperType::BasisWithOperatorsType;
	using FieldType              = typename VectorType::value_type;
	using RealType               = typename BasisWithOperatorsType::RealType;
	using ThisType               = MultiPointCorrelations<CorrelationsSkeletonType>;
	using SparseMatrixType       = typename CorrelationsSkeletonType::SparseMatrixType;

public:

	using MatrixType = typename ObserverHelperType::MatrixType;

	MultiPointCorrelations(const CorrelationsSkeletonType& skeleton)
	    : skeleton_(skeleton)
	{ }

	template <typename VectorLikeType>
	typename PsimagLite::EnableIf<PsimagLite::IsVectorLike<VectorLikeType>::True, void>::Type
	operator()(VectorLikeType&         result,
	           const SparseMatrixType& O,
	           SizeType                rows,
	           SizeType                cols,
	           PsimagLite::String      bra,
	           PsimagLite::String      ket)
	{
		assert(rows == cols);
		result.resize(rows);

		SparseMatrixType Og;

		SparseMatrixType identity(O.rows(), O.cols());
		identity.makeDiagonal(O.rows(), 1.0);

		const size_t rowsOver2 = static_cast<size_t>(rows / 2);

		for (SizeType i = 0; i < rowsOver2; ++i)
			result[i] = calcCorrelation_(Og, i, O, identity, bra, ket);

		for (SizeType i = rowsOver2; i < rows; i++)
			result[i] = 0;
	}

private:

	// from i to i+1
	FieldType calcCorrelation_(SparseMatrixType&       O2gt,
	                           SizeType                i,
	                           const SparseMatrixType& O,
	                           const SparseMatrixType& identity,
	                           PsimagLite::String      bra,
	                           PsimagLite::String      ket)
	{

		if (i >= skeleton_.numberOfSites() - 1)
			throw PsimagLite::RuntimeError("calcCorrelation: i must be < sites-1\n");
		ProgramGlobals::FermionOrBosonEnum fermionicSign
		    = ProgramGlobals::FermionOrBosonEnum::BOSON;

		SizeType         ns = i;
		SparseMatrixType O2g;
		if (i == 0) {
			skeleton_.growDirectly(O2gt, O, i, fermionicSign, ns, true);
			const SizeType ptr
			    = skeleton_.dmrgMultiply(O2g, O2gt, identity, fermionicSign, ns);
			FieldType ret = skeleton_.bracket(O2g, fermionicSign, ptr, bra, ket);
			return ret;
		}

		//			if (i==5) {
		const SizeType ptr = skeleton_.dmrgMultiply(O2g, O2gt, O, fermionicSign, ns - 1);
		//			} else {
		//				skeleton_.dmrgMultiply(O2g,O2gt,identity,fermionicSign,ns-1);
		//			}
		O2gt.clear();
		FieldType ret = skeleton_.bracket(O2g, fermionicSign, ptr, bra, ket);
		skeleton_.helper().transform(O2gt, O2g, ns - 1);
		return ret;
	}

	const CorrelationsSkeletonType& skeleton_;
}; // class MultiPointCorrelations
} // namespace Dmrg

/*@}*/
#endif // MULTI_POINT_CORRELATIONS_H
