// SPDX-FileCopyrightText: Copyright (c) 2012-2017, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 4.]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file KronMatrix.h
 *
 *
 */

#ifndef KRON_MATRIX_HEADER_H
#define KRON_MATRIX_HEADER_H

#include "BatchedGemmInclude.hh"
#include "KronConnections.h"
#include <PsimagLite/Concurrency.h>
#include <PsimagLite/LoadBalancerWeights.h>
#include <PsimagLite/Matrix.h>
#include <PsimagLite/Parallelizer.h>
#include <PsimagLite/ProgressIndicator.h>
#include <PsimagLite/PsimagLite.h>

namespace Dmrg {

template <typename InitKronType> class KronMatrix {

	using SparseMatrixType        = typename InitKronType::SparseMatrixType;
	using ComplexOrRealType       = typename SparseMatrixType::value_type;
	using KronConnectionsType     = KronConnections<InitKronType>;
	using MatrixType              = typename KronConnectionsType::MatrixType;
	using VectorType              = typename KronConnectionsType::VectorType;
	using KronLoggerType          = typename KronConnectionsType::KronLoggerType;
	using ArrayOfMatStructType    = typename InitKronType::ArrayOfMatStructType;
	using GenIjPatchType          = typename InitKronType::GenIjPatchType;
	using MatrixDenseOrSparseType = typename ArrayOfMatStructType::MatrixDenseOrSparseType;
	using VectorSizeType          = typename PsimagLite::Vector<SizeType>::Type;
	using BasisType               = typename GenIjPatchType::BasisType;
	using BatchedGemmType         = BATCHED_GEMM<InitKronType>;

public:

	KronMatrix(InitKronType& initKron, PsimagLite::String name)
	    : initKron_(initKron)
	    , kron_logger_(initKron)
	    , progress_("KronMatrix")
	    , batchedGemm_(initKron)
	{
		PsimagLite::String        str((initKron.loadBalance()) ? "true" : "false");
		PsimagLite::OstringStream msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg << "KronMatrix: " << name << " sizes=" << initKron.size(InitKronType::NEW);
		msg << " " << initKron.size(InitKronType::OLD);
		msg << " loadBalance " << str;
		progress_.printline(msgg, std::cout);
	}

	void matrixVectorProduct(VectorType& vout, const VectorType& vin) const
	{
		initKron_.copyIn(vout, vin);

		if (batchedGemm_.enabled()) {
			VectorType& xout = initKron_.xout();
			VectorType  xoutTmp(xout.size(), 0.0);
			batchedGemm_.matrixVector(xoutTmp, initKron_.yin());
			for (SizeType i = 0; i < xoutTmp.size(); ++i)
				xout[i] += xoutTmp[i];

			initKron_.copyOut(vout);
			return;
		}

		KronConnectionsType kc(initKron_, kron_logger_);
		SizeType            threads = PsimagLite::Concurrency::codeSectionParams.npthreads;
		PsimagLite::CodeSectionParams codeSectionParams(threads);

		if (initKron_.loadBalance()) {
			PsimagLite::Parallelizer<KronConnectionsType,
			                         PsimagLite::LoadBalancerWeights>
			    parallelConnections(codeSectionParams);
			parallelConnections.loopCreate(kc, initKron_.weightsOfPatchesNew());
		} else {
			PsimagLite::Parallelizer<KronConnectionsType> parallelConnections(
			    codeSectionParams);
			parallelConnections.loopCreate(kc);
		}

		kc.sync();

		initKron_.copyOut(vout);
	}

private:

	KronMatrix(const KronMatrix&);

	const KronMatrix& operator=(const KronMatrix&);

	InitKronType&                 initKron_;
	mutable KronLoggerType        kron_logger_;
	PsimagLite::ProgressIndicator progress_;
	BatchedGemmType               batchedGemm_;
}; // class KronMatrix

} // namespace PsimagLite

/*@}*/

#endif // KRON_MATRIX_HEADER_H
