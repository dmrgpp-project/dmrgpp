// SPDX-FileCopyrightText: Copyright (c) 2009-2015, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/
/** \file ParallelWftSu2.h
 */

#ifndef DMRG_PARALLEL_WFT_SU2_H
#define DMRG_PARALLEL_WFT_SU2_H

#include "ProgramGlobals.h"
#include <PsimagLite/Concurrency.h>
#include <PsimagLite/PackIndices.h>
#include <PsimagLite/Vector.h>

namespace Dmrg {

template <typename VectorWithOffsetType, typename DmrgWaveStructType, typename LeftRightSuperType>
class ParallelWftSu2 {

	using PackIndicesType            = PsimagLite::PackIndices;
	using ConcurrencyType            = PsimagLite::Concurrency;
	using VectorVectorWithOffsetType = typename PsimagLite::Vector<VectorWithOffsetType>::Type;
	using VectorSizeType             = PsimagLite::Vector<SizeType>::Type;
	using BasisWithOperatorsType     = typename DmrgWaveStructType::BasisWithOperatorsType;
	using SparseMatrixType           = typename BasisWithOperatorsType::SparseMatrixType;
	using SparseElementType          = typename SparseMatrixType::value_type;
	using FactorsType                = typename BasisWithOperatorsType::FactorsType;

public:

	using VectorElementType = typename VectorWithOffsetType::value_type;
	using RealType          = typename PsimagLite::Real<VectorElementType>::Type;

	ParallelWftSu2(VectorWithOffsetType&         psiDest,
	               const VectorWithOffsetType&   psiSrc,
	               const LeftRightSuperType&     lrs,
	               SizeType                      i0,
	               const VectorSizeType&         nk,
	               const DmrgWaveStructType&     dmrgWaveStruct,
	               ProgramGlobals::DirectionEnum dir)
	    : psiDest_(psiDest)
	    , psiSrc_(psiSrc)
	    , lrs_(lrs)
	    , i0_(i0)
	    , nk_(nk)
	    , dmrgWaveStruct_(dmrgWaveStruct)
	    , dir_(dir)
	    , pack1_(0)
	    , pack2_(0)
	{
		dmrgWaveStruct_.getTransform(ProgramGlobals::SysOrEnvEnum::ENVIRON).toSparse(we_);
		dmrgWaveStruct_.getTransform(ProgramGlobals::SysOrEnvEnum::SYSTEM).toSparse(ws_);
		transposeConjugate(wsT_, ws_);
		transposeConjugate(weT_, we_);

		if (dir_ == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) {
			assert(
			    dmrgWaveStruct_.lrs().right().permutationInverse().size()
			    == dmrgWaveStruct_.getTransform(ProgramGlobals::SysOrEnvEnum::ENVIRON)
			           .rows());
			assert(lrs_.left().permutationInverse().size() / volumeOf(nk)
			       == dmrgWaveStruct_.getTransform(ProgramGlobals::SysOrEnvEnum::SYSTEM)
			              .cols());
			pack1_ = new PackIndicesType(lrs.left().permutationInverse().size());
			pack2_ = new PackIndicesType(lrs.left().permutationInverse().size()
			                             / volumeOf(nk));
		} else {
			assert(dmrgWaveStruct_.lrs().left().permutationInverse().size()
			       == dmrgWaveStruct_.getTransform(ProgramGlobals::SysOrEnvEnum::SYSTEM)
			              .rows());
			assert(
			    lrs_.right().permutationInverse().size() / volumeOf(nk)
			    == dmrgWaveStruct_.getTransform(ProgramGlobals::SysOrEnvEnum::ENVIRON)
			           .cols());
			pack1_ = new PackIndicesType(lrs.super().permutationInverse().size()
			                             / lrs.right().permutationInverse().size());
			pack2_ = new PackIndicesType(volumeOf(nk));
		}

		const FactorsType* fptrSE = lrs_.super().getFactors();
		assert(fptrSE);
		const FactorsType& factorsSE = *fptrSE;

		const FactorsType* fptrS = lrs_.left().getFactors();
		assert(fptrS);
		const FactorsType& factorsS = *fptrS;

		const FactorsType* fptrE = lrs_.right().getFactors();
		assert(fptrE);
		const FactorsType& factorsE = *fptrE;

		transposeConjugate(factorsInvSE_, factorsSE);
		transposeConjugate(factorsInvS_, factorsS);
		transposeConjugate(factorsInvE_, factorsE);
	}

	~ParallelWftSu2()
	{
		delete pack1_;
		delete pack2_;
	}

	static SizeType volumeOf(const VectorSizeType& v)
	{
		assert(v.size() > 0);
		SizeType ret = v[0];
		for (SizeType i = 1; i < v.size(); i++)
			ret *= v[i];
		return ret;
	}

	SizeType tasks() const { return psiDest_.effectiveSize(i0_); }

	void doTask(SizeType taskNumber, SizeType)
	{
		SizeType start                       = psiDest_.offset(i0_);
		psiDest_.fastAccess(i0_, taskNumber) = 0.0;
		SizeType xx                          = taskNumber + start;
		assert(dir_ == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM
		       || dir_ == ProgramGlobals::DirectionEnum::EXPAND_ENVIRON);

		if (dir_ == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) {
			SizeType ip    = 0;
			SizeType alpha = 0;
			SizeType kp    = 0;
			SizeType jp    = 0;
			for (int kI = factorsInvSE_.getRowPtr(xx);
			     kI < factorsInvSE_.getRowPtr(xx + 1);
			     kI++) {
				pack1_->unpack(
				    alpha, jp, static_cast<SizeType>(factorsInvSE_.getCol(kI)));
				for (int k2I = factorsInvS_.getRowPtr(alpha);
				     k2I < factorsInvS_.getRowPtr(alpha + 1);
				     k2I++) {
					pack2_->unpack(
					    ip,
					    kp,
					    static_cast<SizeType>(factorsInvS_.getCol(k2I)));
					psiDest_.fastAccess(i0_, taskNumber)
					    += factorsInvSE_.getValue(kI)
					    * factorsInvS_.getValue(k2I)
					    * createAux2b(psiSrc_, ip, kp, jp, wsT_, we_, nk_);
				}
			}
		} else {
			SizeType ip   = 0;
			SizeType beta = 0;
			SizeType kp   = 0;
			SizeType jp   = 0;
			for (int kI = factorsInvSE_.getRowPtr(xx);
			     kI < factorsInvSE_.getRowPtr(xx + 1);
			     kI++) {
				pack1_->unpack(
				    ip, beta, static_cast<SizeType>(factorsInvSE_.getCol(kI)));
				for (int k2I = factorsInvE_.getRowPtr(beta);
				     k2I < factorsInvE_.getRowPtr(beta + 1);
				     k2I++) {
					pack2_->unpack(
					    kp,
					    jp,
					    static_cast<SizeType>(factorsInvE_.getCol(k2I)));
					psiDest_.fastAccess(i0_, taskNumber)
					    += factorsInvSE_.getValue(kI)
					    * factorsInvE_.getValue(k2I)
					    * createAux1b(psiSrc_, ip, kp, jp, ws_, weT_, nk_);
				}
			}
		}
	}

private:

	// This class has pointers, disallow copy ctor and assignment
	template <typename T1, typename T2, typename T3>
	ParallelWftSu2(const ParallelWftSu2<T1, T2, T3>&);

	template <typename T1, typename T2, typename T3>
	ParallelWftSu2& operator=(const ParallelWftSu2<T1, T2, T3>&);

	template <typename SomeVectorType>
	SparseElementType createAux2b(const SomeVectorType&   psiSrc,
	                              SizeType                ip,
	                              SizeType                kp,
	                              SizeType                jp,
	                              const SparseMatrixType& wsT,
	                              const SparseMatrixType& we,
	                              const VectorSizeType&   nk) const
	{
		SizeType nalpha = dmrgWaveStruct_.lrs().left().permutationInverse().size();
		assert(nalpha == wsT.cols());

		const FactorsType* fptrE = dmrgWaveStruct_.lrs().right().getFactors();
		assert(fptrE);
		const FactorsType& factorsE = *fptrE;

		const FactorsType* fptrSE = dmrgWaveStruct_.lrs().super().getFactors();
		assert(fptrSE);
		const FactorsType& factorsSE = *fptrSE;

		SizeType          volumeOfNk = this->volumeOf(nk);
		SparseElementType sum        = 0;

		SizeType kpjp = kp + jp * volumeOfNk;
		assert(kpjp < dmrgWaveStruct_.lrs().right().permutationInverse().size());
		SizeType kpjpx = dmrgWaveStruct_.lrs().right().permutationInverse(kpjp);

		for (int k2I = factorsE.getRowPtr(kpjpx); k2I < factorsE.getRowPtr(kpjpx + 1);
		     k2I++) {
			SizeType beta = factorsE.getCol(k2I);
			for (int k = wsT.getRowPtr(ip); k < wsT.getRowPtr(ip + 1); k++) {
				SizeType alpha = wsT.getCol(k);
				for (int k2 = we.getRowPtr(beta); k2 < we.getRowPtr(beta + 1);
				     k2++) {
					SizeType j = we.getCol(k2);
					SizeType r = alpha + j * nalpha;
					for (int kI = factorsSE.getRowPtr(r);
					     kI < factorsSE.getRowPtr(r + 1);
					     kI++) {
						SizeType x = factorsSE.getCol(kI);
						sum += wsT.getValue(k) * we.getValue(k2)
						    * psiSrc.slowAccess(x) * factorsSE.getValue(kI)
						    * factorsE.getValue(k2I);
					}
				}
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
	                              const SparseMatrixType& weT,
	                              const VectorSizeType&   nk) const
	{
		SizeType volumeOfNk = volumeOf(nk);
		SizeType ni
		    = dmrgWaveStruct_.getTransform(ProgramGlobals::SysOrEnvEnum::SYSTEM).cols();
		SizeType nip
		    = dmrgWaveStruct_.lrs().left().permutationInverse().size() / volumeOfNk;

		const FactorsType* fptrS = dmrgWaveStruct_.lrs().left().getFactors();
		assert(fptrS);
		const FactorsType& factorsS = *fptrS;

		const FactorsType* fptrSE = dmrgWaveStruct_.lrs().super().getFactors();
		assert(fptrSE);
		const FactorsType& factorsSE = *fptrSE;

		SparseElementType sum = 0;

		SizeType ipkp = ip + kp * nip;
		for (int k2I = factorsS.getRowPtr(ipkp); k2I < factorsS.getRowPtr(ipkp + 1);
		     k2I++) {
			SizeType alpha = factorsS.getCol(k2I);
			for (int k = ws.getRowPtr(alpha); k < ws.getRowPtr(alpha + 1); k++) {
				SizeType i = ws.getCol(k);
				for (int k2 = weT.getRowPtr(jp); k2 < weT.getRowPtr(jp + 1); k2++) {
					SizeType j = weT.getCol(k2);
					SizeType r = i + j * ni;
					for (int kI = factorsSE.getRowPtr(r);
					     kI < factorsSE.getRowPtr(r + 1);
					     kI++) {
						SizeType x = factorsSE.getCol(kI);
						sum += ws.getValue(k) * weT.getValue(k2)
						    * psiSrc.slowAccess(x) * factorsSE.getValue(kI)
						    * factorsS.getValue(k2I);
					}
				}
			}
		}

		return sum;
	}

	VectorWithOffsetType&         psiDest_;
	const VectorWithOffsetType&   psiSrc_;
	const LeftRightSuperType&     lrs_;
	SizeType                      i0_;
	const VectorSizeType&         nk_;
	const DmrgWaveStructType&     dmrgWaveStruct_;
	ProgramGlobals::DirectionEnum dir_;
	SparseMatrixType              ws_;
	SparseMatrixType              we_;
	PackIndicesType*              pack1_;
	PackIndicesType*              pack2_;
	SparseMatrixType              wsT_;
	SparseMatrixType              weT_;
	FactorsType                   factorsInvSE_;
	FactorsType                   factorsInvS_;
	FactorsType                   factorsInvE_;
}; // class ParallelWftSu2
} // namespace Dmrg

/*@}*/
#endif // DMRG_PARALLEL_WFT_SU2_H
