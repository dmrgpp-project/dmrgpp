
// SPDX-FileCopyrightText: Copyright (c) 2009-2015, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [DMRG++ authors; see AUTHORS.]

#ifndef DENSITY_MATRIX_LOCAL_H
#define DENSITY_MATRIX_LOCAL_H
#include "BlockDiagonalMatrix.h"
#include "DensityMatrixBase.h"
#include "DiagBlockDiagMatrix.h"
#include "ParallelDensityMatrix.h"
#include <PsimagLite/Concurrency.h>
#include <PsimagLite/NoPthreads.h>
#include <PsimagLite/Parallelizer.h>
#include <PsimagLite/ProgressIndicator.h>
#include <PsimagLite/TypeToString.h>

namespace Dmrg {

template <typename TargetingType>
class DensityMatrixLocal : public DensityMatrixBase<TargetingType> {

	using BaseType                 = DensityMatrixBase<TargetingType>;
	using LeftRightSuperType       = typename TargetingType::LeftRightSuperType;
	using BasisWithOperatorsType   = typename TargetingType::BasisWithOperatorsType;
	using BasisType                = typename BasisWithOperatorsType::BasisType;
	using SparseMatrixType         = typename BasisWithOperatorsType::SparseMatrixType;
	using TargetVectorType         = typename TargetingType::VectorWithOffsetType;
	using DensityMatrixElementType = typename TargetingType::TargetVectorType::value_type;
	using ConcurrencyType          = PsimagLite::Concurrency;
	using ProgressIndicatorType    = PsimagLite::ProgressIndicator;
	using RealType                 = typename PsimagLite::Real<DensityMatrixElementType>::Type;
	using ParamsType               = typename DensityMatrixBase<TargetingType>::Params;

public:

	using BlockDiagonalMatrixType   = typename BaseType::BlockDiagonalMatrixType;
	using BuildingBlockType         = typename BlockDiagonalMatrixType::BuildingBlockType;
	using ParallelDensityMatrixType = ParallelDensityMatrix<BlockDiagonalMatrixType,
	                                                        BasisWithOperatorsType,
	                                                        TargetVectorType>;
	using ParallelizerType          = PsimagLite::Parallelizer<ParallelDensityMatrixType>;
	using VectorVectorVectorWithOffsetType =
	    typename TargetingType::VectorVectorVectorWithOffsetType;

	DensityMatrixLocal(const TargetingType&      target,
	                   const LeftRightSuperType& lrs,
	                   const ParamsType&         p)
	    : progress_("DensityMatrixLocal")
	    , data_((p.direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) ? lrs.left()
	                                                                          : lrs.right())
	    , direction_(p.direction)
	    , debug_(p.debug)
	{
		{
			PsimagLite::OstringStream                     msgg(std::cout.precision());
			PsimagLite::OstringStream::OstringStreamType& msg = msgg();
			msg << "Init partition for all targets";
			progress_.printline(msgg, std::cout);
		}

		const BasisWithOperatorsType& pBasis
		    = (p.direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) ? lrs.left()
		                                                                    : lrs.right();

		const BasisWithOperatorsType& pBasisSummed
		    = (p.direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) ? lrs.right()
		                                                                    : lrs.left();

		// loop over all partitions:
		for (SizeType m = 0; m < pBasis.partition() - 1; m++) {
			// size of this partition
			SizeType bs = pBasis.partition(m + 1) - pBasis.partition(m);

			// density matrix block for this partition:
			BuildingBlockType matrixBlock(bs, bs);

			// weight of the ground state:
			RealType w = target.gsWeight();

			// if we are to target the ground state do it now:
			if (target.includeGroundStage()) {
				const VectorVectorVectorWithOffsetType& psi = target.psiConst();
				const SizeType                          nsectors = psi.size();

				for (SizeType sectorIndex = 0; sectorIndex < nsectors;
				     ++sectorIndex) {
					const SizeType nexcited = psi[sectorIndex].size();

					for (SizeType excitedIndex = 0; excitedIndex < nexcited;
					     ++excitedIndex) {

						initPartition(matrixBlock,
						              pBasis,
						              m,
						              *(psi[sectorIndex][excitedIndex]),
						              pBasisSummed,
						              lrs.super(),
						              p.direction,
						              w);
					}
				}
			}

			// target all other states if any:
			for (SizeType ix = 0; ix < target.size(); ++ix) {
				RealType wnorm = target.normSquared(ix);
				if (fabs(wnorm) < 1e-6)
					continue;
				RealType w = target.weight(ix) / wnorm;
				initPartition(matrixBlock,
				              pBasis,
				              m,
				              target(ix),
				              pBasisSummed,
				              lrs.super(),
				              p.direction,
				              w);
			}

			// set this matrix block into data_
			data_.setBlock(m, pBasis.partition(m), matrixBlock);
		}

		{
			PsimagLite::OstringStream                     msgg(std::cout.precision());
			PsimagLite::OstringStream::OstringStreamType& msg = msgg();
			msg << "Done with init partition";
			progress_.printline(msgg, std::cout);
		}
	}

	const BlockDiagonalMatrixType& operator()() override { return data_; }

	void diag(typename PsimagLite::Vector<RealType>::Type& eigs, char jobz) override
	{
		DiagBlockDiagMatrix<BlockDiagonalMatrixType>::diagonalise(data_, eigs, jobz);
	}

	friend std::ostream& operator<<(std::ostream& os, const DensityMatrixLocal& dm)
	{
		for (SizeType m = 0; m < dm.data_.blocks(); ++m) {
			SizeType ne = dm.pBasis_.electrons(dm.pBasis_.partition(m));
			os << " ne=" << ne << "\n";
			os << dm.data_(m) << "\n";
		}

		return os;
	}

private:

	void initPartition(BuildingBlockType&            matrixBlock,
	                   BasisWithOperatorsType const& pBasis,
	                   SizeType                      m,
	                   const TargetVectorType&       v,
	                   BasisWithOperatorsType const& pBasisSummed,
	                   BasisType const&              pSE,
	                   ProgramGlobals::DirectionEnum direction,
	                   RealType                      weight)
	{
		ParallelDensityMatrixType helperDm(
		    v, pBasis, pBasisSummed, pSE, direction, m, weight, matrixBlock);
		ParallelizerType threadedDm(ConcurrencyType::codeSectionParams);
		threadedDm.loopCreate(helperDm);
	}

	ProgressIndicatorType         progress_;
	BlockDiagonalMatrixType       data_;
	ProgramGlobals::DirectionEnum direction_;
	bool                          debug_;
}; // class DensityMatrixLocal

} // namespace Dmrg

#endif
