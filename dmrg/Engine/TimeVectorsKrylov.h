// SPDX-FileCopyrightText: Copyright (c) 2009-2012, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 2.0.0]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file TimeVectorsKrylov.h
 *
 *
 */

#ifndef TIME_VECTORS_KRYLOV
#define TIME_VECTORS_KRYLOV
#include "KrylovHelper.h"
#include "ParallelTriDiag.h"
#include "TimeVectorsBase.h"
#include <PsimagLite/NoPthreadsNg.h>
#include <PsimagLite/Parallelizer.h>
#include <iostream>
#include <vector>

namespace Dmrg {

template <typename TargetParamsType,
          typename ModelType,
          typename WaveFunctionTransfType,
          typename LanczosSolverType,
          typename VectorWithOffsetType_>
class TimeVectorsKrylov : public TimeVectorsBase<TargetParamsType,
                                                 ModelType,
                                                 WaveFunctionTransfType,
                                                 LanczosSolverType,
                                                 VectorWithOffsetType_> {

public:

	using VectorWithOffsetType   = VectorWithOffsetType_;
	using BaseType               = TimeVectorsBase<TargetParamsType,
	                                               ModelType,
	                                               WaveFunctionTransfType,
	                                               LanczosSolverType,
	                                               VectorWithOffsetType>;
	using ThisType               = TimeVectorsKrylov<TargetParamsType,
	                                                 ModelType,
	                                                 WaveFunctionTransfType,
	                                                 LanczosSolverType,
	                                                 VectorWithOffsetType>;
	using PairType               = typename BaseType::PairType;
	using RealType               = typename TargetParamsType::RealType;
	using VectorRealType         = typename PsimagLite::Vector<RealType>::Type;
	using ModelHelperType        = typename ModelType::ModelHelperType;
	using LeftRightSuperType     = typename ModelHelperType::LeftRightSuperType;
	using BasisWithOperatorsType = typename LeftRightSuperType::BasisWithOperatorsType;
	using SparseMatrixType       = typename BasisWithOperatorsType::SparseMatrixType;
	using ComplexOrRealType      = typename SparseMatrixType::value_type;
	using ParallelTriDiagType
	    = ParallelTriDiag<ModelType, LanczosSolverType, VectorWithOffsetType>;
	using MatrixComplexOrRealType    = typename ParallelTriDiagType::MatrixComplexOrRealType;
	using VectorType                 = typename ParallelTriDiagType::TargetVectorType;
	using VectorMatrixFieldType      = typename ParallelTriDiagType::VectorMatrixFieldType;
	using TridiagonalMatrixType      = typename LanczosSolverType::TridiagonalMatrixType;
	using InputValidatorType         = typename ModelType::InputValidatorType;
	using VectorVectorWithOffsetType = typename PsimagLite::Vector<VectorWithOffsetType*>::Type;
	using VectorVectorRealType       = typename PsimagLite::Vector<VectorRealType>::Type;

	struct Action {

		using VectorRealType = typename ThisType::VectorRealType;
	};

	struct TypeWrapper {

		using MatrixComplexOrRealType = typename ThisType::MatrixComplexOrRealType;
		using VectorWithOffsetType    = typename ThisType::VectorWithOffsetType;
		using SolverParamsType        = typename ModelType::SolverParamsType;
	};

	using KrylovHelperType = KrylovHelper<Action, TypeWrapper>;

	TimeVectorsKrylov(const TargetParamsType&       tstStruct,
	                  VectorVectorWithOffsetType&   targetVectors,
	                  const ModelType&              model,
	                  const WaveFunctionTransfType& wft,
	                  const LeftRightSuperType&     lrs,
	                  InputValidatorType&           ioIn)
	    : BaseType(model, lrs, wft, "krylov")
	    , tstStruct_(tstStruct)
	    , targetVectors_(targetVectors)
	    , model_(model)
	    , wft_(wft)
	    , lrs_(lrs)
	    , ioIn_(ioIn)
	    , timeHasAdvanced_(false)
	    , krylovHelper_(model.params(), 0)
	{ }

	void calcTimeVectors(const PsimagLite::Vector<SizeType>::Type& indices,
	                     RealType                                  Eg,
	                     const VectorWithOffsetType&               phi,
	                     const typename BaseType::ExtraData&       extra) override
	{
		const VectorRealType& times = tstStruct_.times();

		if (indices.size() < 2)
			err("TimeVectorsKrylov: indices.size() must be greater than 1\n");

		// Only for RIXS Dynamics wft and (if needed) advance
		if (extra.wftAndAdvanceIfNeeded) {
			SizeType advance
			    = (timeHasAdvanced_) ? indices[indices.size() - 1] : indices[0];
			VectorWithOffsetType phiNew;
			assert(targetVectors_[advance]);
			if (targetVectors_[advance]->size() > 0) {
				SizeType numberOfSites = lrs_.super().block().size();
				if (extra.block[0] != 0 && extra.block[0] != numberOfSites - 1) {
					BaseType::wftHelper().wftOneVector(
					    phiNew, *targetVectors_[advance], extra.block[0]);

					*targetVectors_[indices[0]] = phiNew;
				} else {
					if (timeHasAdvanced_) {
						assert(indices[0] != advance);
						*targetVectors_[indices[0]]
						    = *targetVectors_[advance];
					}
				}
			}
		}

		const SizeType n = indices.size();
		if (this->currentTimeStep() == 0 && tstStruct_.noOperator()
		    && tstStruct_.skipTimeZero()) {
			for (SizeType i = 0; i < n; ++i) {
				const SizeType ii   = indices[i];
				*targetVectors_[ii] = phi;
			}
		}

		const VectorWithOffsetType* ptr0 = targetVectors_[indices[0]];
		const VectorWithOffsetType* ptr1 = &phi;
		if (ptr0 != ptr1)
			*targetVectors_[indices[0]] = phi;

		if (times.size() == 1 && fabs(times[0]) < 1e-10)
			return;

		VectorMatrixFieldType V(phi.sectors());
		VectorMatrixFieldType T(phi.sectors());

		typename PsimagLite::Vector<SizeType>::Type steps(phi.sectors());

		triDiag(phi, T, V, steps);

		VectorVectorRealType eigs(phi.sectors());

		for (SizeType ii = 0; ii < phi.sectors(); ii++)
			PsimagLite::diag(T[ii], eigs[ii], 'V');

		calcTargetVectors(indices, phi, T, V, Eg, eigs, steps);

		// checkNorms();
		if (extra.isLastCall)
			timeHasAdvanced_ = false;
	}

	void timeHasAdvanced() override
	{
		timeHasAdvanced_ = true;
		this->advanceCurrentTime(tstStruct_.tau());
	}

private:

	//! Do not normalize states here, it leads to wrong results (!)
	void calcTargetVectors(typename PsimagLite::Vector<SizeType>::Type indices,
	                       const VectorWithOffsetType&                 phi,
	                       const VectorMatrixFieldType&                T,
	                       const VectorMatrixFieldType&                V,
	                       RealType                                    Eg,
	                       const VectorVectorRealType&                 eigs,
	                       typename PsimagLite::Vector<SizeType>::Type steps)
	{
		const VectorRealType& times = tstStruct_.times();

		for (SizeType i = 1; i < indices.size(); ++i) {
			const SizeType ii = indices[i];
			assert(ii < targetVectors_.size());
			*targetVectors_[ii] = phi;
			// Only time differences here (i.e. extra.times[i] not
			// extra.times[i]+currentTime_)
			calcTargetVector(*targetVectors_[ii], phi, T, V, Eg, eigs, steps, i, times);
		}
	}

	void calcTargetVector(VectorWithOffsetType&                       v,
	                      const VectorWithOffsetType&                 phi,
	                      const VectorMatrixFieldType&                T,
	                      const VectorMatrixFieldType&                V,
	                      RealType                                    Eg,
	                      const VectorVectorRealType&                 eigs,
	                      typename PsimagLite::Vector<SizeType>::Type steps,
	                      SizeType                                    timeIndex,
	                      const VectorRealType&                       times)
	{
		v = phi;
		for (SizeType ii = 0; ii < phi.sectors(); ++ii) {
			const RealType        time          = times[timeIndex];
			const RealType        timeDirection = tstStruct_.timeDirection();
			const VectorRealType& eigsii        = eigs[ii];
			auto                  action = [eigsii, Eg, time, timeDirection](SizeType k)
			{
				RealType          tmp = (eigsii[k] - Eg) * time * timeDirection;
				ComplexOrRealType c   = 0.0;
				PsimagLite::expComplexOrReal(c, -tmp);
				return c;
			};

			SizeType   i0 = phi.sector(ii);
			VectorType r;
			calcTargetVector(r, phi, T[ii], V[ii], action, steps[ii], i0);
			v.setDataInSector(r, i0);
		}
	}

	template <typename SomeLambdaType>
	void calcTargetVector(VectorType&                    r,
	                      const VectorWithOffsetType&    phi,
	                      const MatrixComplexOrRealType& T,
	                      const MatrixComplexOrRealType& V,
	                      const SomeLambdaType&          action,
	                      SizeType                       steps,
	                      SizeType                       i0)
	{
		SizeType n2 = steps;
		SizeType n  = V.rows();
		if (T.cols() != T.rows())
			throw PsimagLite::RuntimeError("T is not square\n");
		if (V.cols() != T.cols())
			throw PsimagLite::RuntimeError("V is not nxn2\n");
		// for (SizeType j=0;j<v.size();j++) v[j] = 0; <-- harmful if v is sparse
		ComplexOrRealType zone  = 1.0;
		ComplexOrRealType zzero = 0.0;

		// check1(phi,i0);
		// check2(T,V,phi,n2,i0);
		VectorType tmp(n2);
		r.resize(n2);
		krylovHelper_.calcR(r, action, T, V, phi, steps, i0);
		psimag::BLAS::GEMV(
		    'N', n2, n2, zone, &(T(0, 0)), n2, &(r[0]), 1, zzero, &(tmp[0]), 1);
		r.resize(n);
		psimag::BLAS::GEMV(
		    'N', n, n2, zone, &(V(0, 0)), n, &(tmp[0]), 1, zzero, &(r[0]), 1);
	}

	void triDiag(const VectorWithOffsetType&                  phi,
	             VectorMatrixFieldType&                       T,
	             VectorMatrixFieldType&                       V,
	             typename PsimagLite::Vector<SizeType>::Type& steps)
	{
		using ParallelizerType = PsimagLite::NoPthreadsNg<ParallelTriDiagType>;
		ParallelizerType threadedTriDiag(PsimagLite::CodeSectionParams(1));

		ParallelTriDiagType helperTriDiag(
		    phi, T, V, steps, lrs_, this->time(), model_, ioIn_);

		threadedTriDiag.loopCreate(helperTriDiag);
	}

	const TargetParamsType&       tstStruct_;
	VectorVectorWithOffsetType&   targetVectors_;
	const ModelType&              model_;
	const WaveFunctionTransfType& wft_;
	const LeftRightSuperType&     lrs_;
	InputValidatorType&           ioIn_;
	bool                          timeHasAdvanced_;
	KrylovHelperType              krylovHelper_;
}; // class TimeVectorsKrylov
} // namespace Dmrg
/*@}*/
#endif
