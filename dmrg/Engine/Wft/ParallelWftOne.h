// SPDX-FileCopyrightText: Copyright (c) 2009-2015, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/
/** \file ParallelWftOne.h
 */

#ifndef DMRG_PARALLEL_WFT_ONE_H
#define DMRG_PARALLEL_WFT_ONE_H

#include "OffsetVector.hpp"
#include "ProgramGlobals.h"
#include <PsimagLite/Concurrency.h>
#include <PsimagLite/PackIndices.h>
#include <PsimagLite/Vector.h>

namespace Dmrg {

template <typename VectorWithOffsetType, typename DmrgWaveStructType, typename OneSiteSpacesType>
class ParallelWftOne {

	using PackIndicesType            = PsimagLite::PackIndices;
	using ConcurrencyType            = PsimagLite::Concurrency;
	using VectorVectorWithOffsetType = typename PsimagLite::Vector<VectorWithOffsetType>::Type;
	using VectorSizeType             = typename DmrgWaveStructType::VectorSizeType;
	using BasisWithOperatorsType     = typename DmrgWaveStructType::BasisWithOperatorsType;
	using SparseMatrixType           = typename BasisWithOperatorsType::SparseMatrixType;
	using SparseElementType          = typename SparseMatrixType::value_type;
	using LeftRightSuperType         = typename DmrgWaveStructType::LeftRightSuperType;
	using OffsetVectorType           = OffsetVector<SparseElementType>;
	using OffsetVectorBaseType       = OffsetVectorBase<SparseElementType>;

public:

	using VectorElementType = typename VectorWithOffsetType::value_type;
	using RealType          = typename PsimagLite::Real<VectorElementType>::Type;

	ParallelWftOne(VectorWithOffsetType&       psiDest,
	               const VectorWithOffsetType& psiSrc,
	               const LeftRightSuperType&   lrs,
	               SizeType                    i0,
	               const OneSiteSpacesType&    oneSiteSpaces,
	               const DmrgWaveStructType&   dmrgWaveStruct)
	    : psiDest_(psiDest)
	    , psiSrc_(psiSrc)
	    , lrs_(lrs)
	    , i0_(i0)
	    , oneSiteSpaces_(oneSiteSpaces)
	    , dmrgWaveStruct_(dmrgWaveStruct)
	    , pack1_(0)
	    , pack2_(0)
	{
		dmrgWaveStruct_.getTransform(ProgramGlobals::SysOrEnvEnum::ENVIRON).toSparse(we_);
		dmrgWaveStruct_.getTransform(ProgramGlobals::SysOrEnvEnum::SYSTEM).toSparse(ws_);
		transposeConjugate(wsT_, ws_);
		transposeConjugate(weT_, we_);
		SizeType vOfNk                             = oneSiteSpaces.hilbertMain(); // CHECK!
		typename ProgramGlobals::DirectionEnum dir = oneSiteSpaces.direction();

		if (dir == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) {
			assert(
			    dmrgWaveStruct_.lrs().right().permutationInverse().size()
			    == dmrgWaveStruct_.getTransform(ProgramGlobals::SysOrEnvEnum::ENVIRON)
			           .rows());
			assert(lrs_.left().permutationInverse().size() / vOfNk
			       == dmrgWaveStruct_.getTransform(ProgramGlobals::SysOrEnvEnum::SYSTEM)
			              .cols());
			pack1_ = new PackIndicesType(lrs.left().permutationInverse().size());
			pack2_
			    = new PackIndicesType(lrs.left().permutationInverse().size() / vOfNk);
		} else {
			assert(dmrgWaveStruct_.lrs().left().permutationInverse().size()
			       == dmrgWaveStruct_.getTransform(ProgramGlobals::SysOrEnvEnum::SYSTEM)
			              .rows());
			assert(
			    lrs_.right().permutationInverse().size() / vOfNk
			    == dmrgWaveStruct_.getTransform(ProgramGlobals::SysOrEnvEnum::ENVIRON)
			           .cols());
			pack1_ = new PackIndicesType(lrs.super().permutationInverse().size()
			                             / lrs.right().permutationInverse().size());
			pack2_ = new PackIndicesType(vOfNk);
		}
	}

	~ParallelWftOne()
	{
		delete pack1_;
		delete pack2_;
	}

	SizeType tasks() const { return psiDest_.effectiveSize(i0_); }

	void doTask(SizeType taskNumber, SizeType)
	{
		SizeType                               start = psiDest_.offset(i0_);
		typename ProgramGlobals::DirectionEnum dir   = oneSiteSpaces_.direction();

		if (dir == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) {
			SizeType ip    = 0;
			SizeType alpha = 0;
			SizeType kp    = 0;
			SizeType jp    = 0;
			pack1_->unpack(
			    alpha, jp, (SizeType)lrs_.super().permutation(taskNumber + start));
			pack2_->unpack(ip, kp, (SizeType)lrs_.left().permutation(alpha));
			psiDest_.fastAccess(i0_, taskNumber)
			    = createAux2b(psiSrc_, ip, kp, jp, wsT_, we_);
		} else {
			SizeType ip   = 0;
			SizeType beta = 0;
			SizeType kp   = 0;
			SizeType jp   = 0;
			pack1_->unpack(
			    ip, beta, (SizeType)lrs_.super().permutation(taskNumber + start));
			pack2_->unpack(kp, jp, (SizeType)lrs_.right().permutation(beta));
			psiDest_.fastAccess(i0_, taskNumber)
			    = createAux1b(psiSrc_, ip, kp, jp, ws_, weT_);
		}
	}

private:

	// This class has pointers, disallow copy ctor and assignment
	ParallelWftOne(const ParallelWftOne&);

	ParallelWftOne& operator=(const ParallelWftOne&);

	template <typename SomeVectorType>
	SparseElementType createAux2b(const SomeVectorType&   psiSrc,
	                              SizeType                ip,
	                              SizeType                kp,
	                              SizeType                jp,
	                              const SparseMatrixType& wsT,
	                              const SparseMatrixType& we) const
	{
		SizeType nalpha = dmrgWaveStruct_.lrs().left().permutationInverse().size();
		assert(nalpha == wsT.cols());

		OffsetVectorType src_proxy;
		const auto&      psiSrc_opt = src_proxy.makeOffsetVector(psiSrc);

		SparseElementType sum        = 0;
		SizeType          volumeOfNk = oneSiteSpaces_.hilbertMain(); // CHECK!
		SizeType          beta
		    = dmrgWaveStruct_.lrs().right().permutationInverse(kp + jp * volumeOfNk);

		for (int k = wsT.getRowPtr(ip); k < wsT.getRowPtr(ip + 1); k++) {
			SizeType alpha  = wsT.getCol(k);
			SizeType begink = we.getRowPtr(beta);
			SizeType endk   = we.getRowPtr(beta + 1);
			for (SizeType k2 = begink; k2 < endk; ++k2) {
				SizeType j = we.getCol(k2);
				SizeType x = dmrgWaveStruct_.lrs().super().permutationInverse(
				    alpha + j * nalpha);
				sum += wsT.getValue(k) * we.getValue(k2) * psiSrc_opt.slowAccess(x);
			}
		}

		return sum;
	}

	template <typename SomeVectorType>
	SparseElementType createAux1b(const SomeVectorType&   psiSrc,
	                              SizeType                ip,
	                              SizeType                kp,
	                              SizeType                jp,
	                              const SparseMatrixType& ws,
	                              const SparseMatrixType& weT) const
	{
		SizeType volumeOfNk = oneSiteSpaces_.hilbertMain(); // CHECK!
		SizeType ni
		    = dmrgWaveStruct_.getTransform(ProgramGlobals::SysOrEnvEnum::SYSTEM).cols();
		SizeType nip
		    = dmrgWaveStruct_.lrs().left().permutationInverse().size() / volumeOfNk;
		SizeType alpha = dmrgWaveStruct_.lrs().left().permutationInverse(ip + kp * nip);

		OffsetVectorType src_proxy;
		const auto&      psiSrc_opt = src_proxy.makeOffsetVector(psiSrc);

		SparseElementType sum = 0;

		for (int k = ws.getRowPtr(alpha); k < ws.getRowPtr(alpha + 1); k++) {
			SizeType i = ws.getCol(k);
			for (int k2 = weT.getRowPtr(jp); k2 < weT.getRowPtr(jp + 1); k2++) {
				SizeType j = weT.getCol(k2);
				SizeType x
				    = dmrgWaveStruct_.lrs().super().permutationInverse(i + j * ni);
				sum += ws.getValue(k) * weT.getValue(k2) * psiSrc_opt.slowAccess(x);
			}
		}

		return sum;
	}

	VectorWithOffsetType&       psiDest_;
	const VectorWithOffsetType& psiSrc_;
	const LeftRightSuperType&   lrs_;
	SizeType                    i0_;
	const OneSiteSpacesType&    oneSiteSpaces_;
	const DmrgWaveStructType&   dmrgWaveStruct_;
	SparseMatrixType            we_;
	SparseMatrixType            ws_;
	PackIndicesType*            pack1_;
	PackIndicesType*            pack2_;
	SparseMatrixType            wsT_;
	SparseMatrixType            weT_;
}; // class ParallelWftOne
} // namespace Dmrg

/*@}*/
#endif // DMRG_PARALLEL_WFT_ONE_H
