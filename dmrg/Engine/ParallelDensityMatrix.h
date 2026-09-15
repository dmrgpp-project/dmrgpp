// SPDX-FileCopyrightText: Copyright (c) 2009,-2012 UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 2.0.0]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/
/** \file ParallelDensityMatrix.h
 */

#ifndef PARALLEL_DENSITY_MATRIX_H
#define PARALLEL_DENSITY_MATRIX_H

#include "ProgramGlobals.h"
#include <PsimagLite/Concurrency.h>

namespace Dmrg {

template <typename BlockMatrixType, typename BasisWithOperatorsType, typename TargetVectorType>
class ParallelDensityMatrix {

	using BuildingBlockType        = typename BlockMatrixType::BuildingBlockType;
	using DensityMatrixElementType = typename TargetVectorType::value_type;
	using BasisType                = typename BasisWithOperatorsType::BasisType;
	using ConcurrencyType          = PsimagLite::Concurrency;

public:

	using RealType = typename PsimagLite::Real<DensityMatrixElementType>::Type;

	ParallelDensityMatrix(const TargetVectorType&       target,
	                      const BasisWithOperatorsType& pBasis,
	                      const BasisWithOperatorsType& pBasisSummed,
	                      const BasisType&              pSE,
	                      ProgramGlobals::DirectionEnum direction,
	                      SizeType                      m,
	                      RealType                      weight,
	                      BuildingBlockType&            matrixBlock)
	    : target_(target)
	    , pBasis_(pBasis)
	    , pBasisSummed_(pBasisSummed)
	    , pSE_(pSE)
	    , direction_(direction)
	    , m_(m)
	    , weight_(weight)
	    , matrixBlock_(matrixBlock)
	    , hasMpi_(PsimagLite::Concurrency::hasMpi())
	{ }

	SizeType tasks() const { return pBasis_.partition(m_ + 1) - pBasis_.partition(m_); }

	void doTask(SizeType taskNumber, SizeType threadNum)
	{
		SizeType start  = pBasis_.partition(m_);
		SizeType length = pBasis_.partition(m_ + 1) - start;

		SizeType ieff = taskNumber + start;
		for (SizeType j = 0; j < length; ++j) {
			matrixBlock_(taskNumber, j)
			    += densityMatrixExpand(direction_, ieff, j + start, target_) * weight_;
		}
	}

private:

	DensityMatrixElementType densityMatrixExpand(ProgramGlobals::DirectionEnum direction,
	                                             SizeType                      alpha1,
	                                             SizeType                      alpha2,
	                                             const TargetVectorType&       v)
	{
		if (direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM)
			return densityMatrixExpandSystem(alpha1, alpha2, v);
		else
			return densityMatrixExpandEnviron(alpha1, alpha2, v);
	}

	DensityMatrixElementType
	densityMatrixExpandEnviron(SizeType alpha1, SizeType alpha2, const TargetVectorType& v)
	{
		SizeType                 ns    = pBasisSummed_.size();
		SizeType                 total = pBasisSummed_.size();
		DensityMatrixElementType sum   = 0;

		SizeType x2 = alpha2 * ns;
		SizeType x1 = alpha1 * ns;
		for (SizeType beta = 0; beta < total; beta++) {
			SizeType ii      = pSE_.permutationInverse(beta + x1);
			int      sector1 = v.index2Sector(ii);
			if (sector1 < 0)
				continue;
			SizeType start1 = v.offset(sector1);

			SizeType jj      = pSE_.permutationInverse(beta + x2);
			int      sector2 = v.index2Sector(jj);
			if (sector2 < 0)
				continue;
			SizeType start2 = v.offset(sector2);

			sum += v.fastAccess(sector1, ii - start1)
			    * PsimagLite::conj(v.fastAccess(sector2, jj - start2));
		}
		return sum;
	}

	DensityMatrixElementType
	densityMatrixExpandSystem(SizeType alpha1, SizeType alpha2, const TargetVectorType& v)
	{
		SizeType                 ne    = pBasisSummed_.size();
		SizeType                 ns    = pSE_.size() / ne;
		SizeType                 total = pBasisSummed_.size();
		DensityMatrixElementType sum   = 0;

		SizeType totalNs = total * ns;

		for (SizeType betaNs = 0; betaNs < totalNs; betaNs += ns) {
			SizeType ii      = pSE_.permutationInverse(alpha1 + betaNs);
			int      sector1 = v.index2Sector(ii);
			if (sector1 < 0)
				continue;
			SizeType start1 = v.offset(sector1);

			SizeType jj      = pSE_.permutationInverse(alpha2 + betaNs);
			int      sector2 = v.index2Sector(jj);
			if (sector2 < 0)
				continue;
			SizeType start2 = v.offset(sector2);

			sum += v.fastAccess(sector1, ii - start1)
			    * PsimagLite::conj(v.fastAccess(sector2, jj - start2));
		}

		return sum;
	}

	const TargetVectorType&       target_;
	const BasisWithOperatorsType& pBasis_;
	const BasisWithOperatorsType& pBasisSummed_;
	const BasisType&              pSE_;
	ProgramGlobals::DirectionEnum direction_;
	SizeType                      m_;
	RealType                      weight_;
	BuildingBlockType&            matrixBlock_;
	bool                          hasMpi_;
}; // class ParallelDensityMatrix
} // namespace Dmrg

/*@}*/
#endif // PARALLEL_DENSITY_MATRIX_H
