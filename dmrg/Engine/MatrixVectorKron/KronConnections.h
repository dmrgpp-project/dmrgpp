// SPDX-FileCopyrightText: Copyright (c) 2012, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 2.0.0]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file KronConnections.h
 *
 *
 */

#ifndef KRON_CONNECTIONS_H
#define KRON_CONNECTIONS_H

#include "../KronUtil/MatrixDenseOrSparse.h"
#include "KronLogger.hh"
#include <PsimagLite/Concurrency.h>
#include <PsimagLite/GemmR.h>
#include <PsimagLite/Matrix.h>

namespace Dmrg {

template <typename InitKronType> class KronConnections {

	using SparseMatrixType        = typename InitKronType::SparseMatrixType;
	using ComplexOrRealType       = typename SparseMatrixType::value_type;
	using ArrayOfMatStructType    = typename InitKronType::ArrayOfMatStructType;
	using GenIjPatchType          = typename InitKronType::GenIjPatchType;
	using ConcurrencyType         = PsimagLite::Concurrency;
	using MatrixDenseOrSparseType = typename ArrayOfMatStructType::MatrixDenseOrSparseType;
	using VectorSizeType          = PsimagLite::Vector<SizeType>::Type;

public:

	using KronLoggerType   = KronLogger<typename InitKronType::ModelType>;
	using MatrixType       = PsimagLite::Matrix<ComplexOrRealType>;
	using VectorType       = typename MatrixDenseOrSparseType::VectorType;
	using VectorVectorType = typename PsimagLite::Vector<VectorType>::Type;
	using RealType         = typename InitKronType::RealType;

	KronConnections(InitKronType& initKron, KronLoggerType& kron_logger)
	    : initKron_(initKron)
	    , kron_logger_(kron_logger)
	    , x_(initKron.xout())
	    , y_(initKron.yin())
	{
		kron_logger_.vector(y_);
	}

	SizeType tasks() const { return initKron_.numberOfPatches(InitKronType::NEW); }

	void doTask(SizeType outPatch, SizeType)
	{
		const bool isComplex = PsimagLite::IsComplexNumber<ComplexOrRealType>::True;

		static const bool                    needsPrinting = false;
		PsimagLite::GemmR<ComplexOrRealType> gemmR(
		    needsPrinting, initKron_.gemmRnb(), initKron_.nthreads2());

		SizeType nC      = initKron_.connections();
		SizeType total   = initKron_.numberOfPatches(InitKronType::OLD);
		SizeType offsetX = initKron_.offsetForPatches(InitKronType::NEW, outPatch);
		assert(offsetX < x_.size());
		kron_logger_.one(outPatch);
		for (SizeType inPatch = 0; inPatch < total; ++inPatch) {
			SizeType offsetY = initKron_.offsetForPatches(InitKronType::OLD, inPatch);
			assert(offsetY < y_.size());
			kron_logger_.two(inPatch);
			for (SizeType ic = 0; ic < nC; ++ic) {
				const ArrayOfMatStructType& xiStruct = initKron_.xc(ic);
				const ArrayOfMatStructType& yiStruct = initKron_.yc(ic);

				const bool performTranspose
				    = (initKron_.useLowerPart() && (outPatch < inPatch));

				const MatrixDenseOrSparseType* Amat = performTranspose
				    ? xiStruct(inPatch, outPatch)
				    : xiStruct(outPatch, inPatch);

				const MatrixDenseOrSparseType* Bmat = performTranspose
				    ? yiStruct(inPatch, outPatch)
				    : yiStruct(outPatch, inPatch);

				if (!Amat || !Bmat)
					continue;

				if (!performTranspose)
					initKron_.checks(*Amat, *Bmat, outPatch, inPatch);

				const char opt = performTranspose ? (isComplex ? 'c' : 't') : 'n';
				kron_logger_.three(outPatch, inPatch, ic);
				kronMult(x_,
				         offsetX,
				         y_,
				         offsetY,
				         opt,
				         opt,
				         *Amat,
				         *Bmat,
				         initKron_.denseFlopDiscount(),
				         gemmR);
			}
		}
	}

	void sync() { kron_logger_.sync(); }

private:

	// disable copy ctor
	KronConnections(const KronConnections&);

	// disable assigment operator
	KronConnections& operator=(const KronConnections&);

	const InitKronType& initKron_;
	KronLoggerType&     kron_logger_;
	VectorType&         x_;
	const VectorType&   y_;
}; // class KronConnections

} // namespace PsimagLite

/*@}*/

#endif // KRON_CONNECTIONS_H
