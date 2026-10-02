// SPDX-FileCopyrightText: Copyright (c) 2008 , UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 1.0.0]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file TwoPointCorrelations.h
 *
 *  A class to perform post-processing calculation of TwoPointCorrelations
 *  <state1 | A_i B_j |state2>
 *
 */
#ifndef TWO_POINT_H
#define TWO_POINT_H
#include "ManyPointAction.h"
#include "Parallel2PointCorrelations.h"
#include "ProgramGlobals.h"
#include "SdhsReinterpret.h"
#include "VectorWithOffsets.h" // for operator*
#include <PsimagLite/Concurrency.h>
#include <PsimagLite/CrsMatrix.h>
#include <PsimagLite/GetBraOrKet.h>
#include <PsimagLite/Parallelizer.h>

namespace Dmrg {

template <typename CorrelationsSkeletonType> class TwoPointCorrelations {

public:

	using ObserverHelperType     = typename CorrelationsSkeletonType::ObserverHelperType;
	using VectorType             = typename ObserverHelperType::VectorType;
	using VectorWithOffsetType   = typename ObserverHelperType::VectorWithOffsetType;
	using BasisWithOperatorsType = typename ObserverHelperType::BasisWithOperatorsType;
	using FieldType              = typename VectorType::value_type;
	using RealType               = typename BasisWithOperatorsType::RealType;
	using ThisType               = TwoPointCorrelations<CorrelationsSkeletonType>;
	using BraketType             = typename CorrelationsSkeletonType::BraketType;
	using SparseMatrixType       = typename CorrelationsSkeletonType::SparseMatrixType;
	using MatrixType             = typename ObserverHelperType::MatrixType;
	using Parallel2PointCorrelationsType = Parallel2PointCorrelations<ThisType>;
	using PairType                       = typename Parallel2PointCorrelationsType::PairType;
	using SdhsReinterpretType            = SdhsReinterpret<BraketType>;
	using ManyPointActionType            = ManyPointAction;

	TwoPointCorrelations(const CorrelationsSkeletonType& skeleton)
	    : skeleton_(skeleton)
	{ }

	void operator()(PsimagLite::Matrix<FieldType>&     w,
	                const BraketType&                  braket,
	                ProgramGlobals::FermionOrBosonEnum fermionicSign,
	                const PsimagLite::GetBraOrKet&     bra,
	                const PsimagLite::GetBraOrKet&     ket,
	                const ManyPointActionType&         action) const
	{
		SizeType rows = w.n_row();
		SizeType cols = w.n_col();

		typename PsimagLite::Vector<PairType>::Type pairs;
		for (SizeType i = 0; i < rows; i++) {
			for (SizeType j = i; j < cols; j++) {
				if (i > j)
					continue;
				if (!action(i, j))
					continue;
				pairs.push_back(PairType(i, j));
			}
		}

		using ParallelizerType = PsimagLite::Parallelizer<Parallel2PointCorrelationsType>;
		ParallelizerType threaded2Points(PsimagLite::Concurrency::codeSectionParams);

		Parallel2PointCorrelationsType helper2Points(
		    w, *this, pairs, braket, fermionicSign, bra, ket);

		threaded2Points.loopCreate(helper2Points);
	}

	// Return the vector: O1 * O2 |psi>
	// where |psi> is the g.s.
	// Note1: O1 is applied to site i and O2 is applied to site j
	// Note2: O1 and O2 operators must commute or anti-commute (set fermionicSign accordingly)
	FieldType calcCorrelation(SizeType                           i,
	                          SizeType                           j,
	                          const BraketType&                  braket,
	                          ProgramGlobals::FermionOrBosonEnum fermionicSign,
	                          const PsimagLite::GetBraOrKet&     bra,
	                          const PsimagLite::GetBraOrKet&     ket) const
	{
		FieldType           c = 0;
		SdhsReinterpretType sdhs(braket, { i, j });

		if (sdhs.forbidden())
			return sdhs.forbiddenValue();

		const SparseMatrixType& O1 = sdhs.op(0).getCRS();
		const SparseMatrixType& O2 = sdhs.op(1).getCRS();

		const RealType fsign
		    = (fermionicSign == ProgramGlobals::FermionOrBosonEnum::BOSON) ? 1 : -1;
		if (i == j) {
			SizeType replacementSite = (i == 0) ? 1 : i - 1;
			SizeType rowsForIdent    = braket.model().hilbertSize(replacementSite);
			c                        = calcDiagonalCorrelation(
                            i, O1, O2, rowsForIdent, fermionicSign, bra, ket);
		} else if (i > j) {
			c = fsign * calcCorrelation_(j, i, O2, O1, fermionicSign, bra, ket);
		} else {
			c = calcCorrelation_(i, j, O1, O2, fermionicSign, bra, ket);
		}

		return c;
	}

private:

	FieldType calcDiagonalCorrelation(SizeType                i,
	                                  const SparseMatrixType& O1,
	                                  const SparseMatrixType& O2,
	                                  SizeType                rowsForIdent,
	                                  ProgramGlobals::FermionOrBosonEnum,
	                                  const PsimagLite::GetBraOrKet& bra,
	                                  const PsimagLite::GetBraOrKet& ket) const
	{
		SparseMatrixType ident = identity(rowsForIdent);

		SparseMatrixType O2new = O1 * O2;
		if (i == 0)
			return calcCorrelation_(0,
			                        1,
			                        O2new,
			                        ident,
			                        ProgramGlobals::FermionOrBosonEnum::BOSON,
			                        bra,
			                        ket);

		return calcCorrelation_(
		    i - 1, i, ident, O2new, ProgramGlobals::FermionOrBosonEnum::BOSON, bra, ket);
	}

	FieldType calcCorrelation_(SizeType                           i,
	                           SizeType                           j,
	                           const SparseMatrixType&            O1,
	                           const SparseMatrixType&            O2,
	                           ProgramGlobals::FermionOrBosonEnum fermionicSign,
	                           const PsimagLite::GetBraOrKet&     bra,
	                           const PsimagLite::GetBraOrKet&     ket) const
	{

		if (i >= j)
			err("Observer::calcCorrelation_(...): i must be smaller than j\n");

		const ObserverHelperType& helper = skeleton_.helper();
		SparseMatrixType          O1m, O2m;
		skeleton_.createWithModification(O1m, O1, 'n');
		skeleton_.createWithModification(O2m, O2, 'n');

		if (j == skeleton_.numberOfSites() - 1) {
			if (i == j - 1) {
				const SizeType ptr = j - 2;
				SizeType ni = helper.leftRightSuper(ptr).left().size() / O1m.rows();

				SparseMatrixType O1g;
				O1g.makeDiagonal(ni, 1.0);

				return skeleton_.bracketRightCorner(
				    O1g, O1m, O2m, fermionicSign, ptr, bra, ket);
			}

			SparseMatrixType O1g;
			skeleton_.growDirectly(O1g, O1m, i, fermionicSign, j - 1, false);
			// j - 2 below is the pointer
			return skeleton_.bracketRightCorner(
			    O1g, O2m, fermionicSign, j - 2, bra, ket);
		}

		SparseMatrixType O1g, O2g;
		SizeType         ns = j - 1;

		skeleton_.growDirectly(O1g, O1m, i, fermionicSign, ns, true);
		const SizeType ptr = skeleton_.dmrgMultiply(O2g, O1g, O2m, fermionicSign, ns);

		return skeleton_.bracket(
		    O2g, ProgramGlobals::FermionOrBosonEnum::BOSON, ptr, bra, ket);
	}

	static SparseMatrixType identity(SizeType n)
	{
		SparseMatrixType ret(n, n);
		ret.makeDiagonal(n, 1.0);
		return ret;
	}

	const CorrelationsSkeletonType& skeleton_;
}; // class TwoPointCorrelations
} // namespace Dmrg

/*@}*/
#endif // TWO_POINT_H
