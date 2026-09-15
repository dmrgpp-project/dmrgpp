// SPDX-FileCopyrightText: Copyright (c) 2009-2015, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/
/** \file ParallelWftMany.h
 */

#ifndef DMRG_PARALLEL_WFT_MANY_H
#define DMRG_PARALLEL_WFT_MANY_H

#include <PsimagLite/Concurrency.h>
#include <PsimagLite/Vector.h>

namespace Dmrg {

template <typename VectorWithOffsetType,
          typename WaveFunctionTransfType,
          typename LeftRightSuperType>
class ParallelWftMany {

	using ConcurrencyType            = PsimagLite::Concurrency;
	using VectorVectorWithOffsetType = typename PsimagLite::Vector<VectorWithOffsetType>::Type;

public:

	using VectorElementType = typename VectorWithOffsetType::value_type;
	using RealType          = typename PsimagLite::Real<VectorElementType>::Type;

	ParallelWftMany(VectorVectorWithOffsetType&   targetVectors,
	                SizeType                      nk,
	                const WaveFunctionTransfType& wft,
	                const LeftRightSuperType&     lrs)
	    : targetVectors_(targetVectors)
	    , nk_(nk)
	    , wft_(wft)
	    , lrs_(lrs)
	{ }

	void thread_function_(SizeType threadNum,
	                      SizeType blockSize,
	                      SizeType total,
	                      ConcurrencyType::MutexType*)
	{
		SizeType nk        = nk_;
		SizeType mpiRank   = PsimagLite::MPI::commRank(PsimagLite::MPI::COMM_WORLD);
		SizeType npthreads = PsimagLite::Concurrency::npthreads;

		ConcurrencyType::mpiDisableIfNeeded(mpiRank, blockSize, "ParallelWftMany", total);

		for (SizeType p = 0; p < blockSize; p++) {
			SizeType ix = (threadNum + npthreads * mpiRank) * blockSize + p + 1;
			if (ix >= targetVectors_.size())
				break;
			VectorWithOffsetType phiNew = targetVectors_[0];
			wft_.setInitialVector(phiNew, targetVectors_[ix], lrs_, nk);
			targetVectors_[ix] = phiNew;
		}
	}

private:

	VectorVectorWithOffsetType&   targetVectors_;
	SizeType                      nk_;
	const WaveFunctionTransfType& wft_;
	const LeftRightSuperType&     lrs_;
}; // class ParallelWftMany
} // namespace Dmrg

/*@}*/
#endif // DMRG_PARALLEL_WFT_MANY_H
