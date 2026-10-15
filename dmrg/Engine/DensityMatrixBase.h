// SPDX-FileCopyrightText: Copyright (c) 2009-2017, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 4.]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file DensityMatrixBase.h
 *
 *
 *
 */
#ifndef DENSITY_MATRIX_BASE_H
#define DENSITY_MATRIX_BASE_H

#include "BlockDiagonalMatrix.h"
#include <PsimagLite/PsimagLite.h>

namespace Dmrg {
template <typename TargetingType> class DensityMatrixBase {

public:

	using BasisWithOperatorsType   = typename TargetingType::BasisWithOperatorsType;
	using BasisType                = typename BasisWithOperatorsType::BasisType;
	using SparseMatrixType         = typename BasisWithOperatorsType::SparseMatrixType;
	using DensityMatrixElementType = typename TargetingType::TargetVectorType::value_type;
	using RealType                 = typename PsimagLite::Real<DensityMatrixElementType>::Type;
	using VectorRealType           = typename PsimagLite::Vector<RealType>::Type;
	using MatrixType               = PsimagLite::Matrix<DensityMatrixElementType>;
	using BlockDiagonalMatrixType  = BlockDiagonalMatrix<MatrixType>;

	struct Params {

		Params(bool                          u,
		       ProgramGlobals::DirectionEnum d,
		       bool                          de,
		       bool                          enablePersistentSvd_,
		       bool                          serialSvd_)
		    : useSvd(u)
		    , direction(d)
		    , debug(de)
		    , enablePersistentSvd(enablePersistentSvd_)
		    , serialSvd(serialSvd_)
		{ }

		bool                          useSvd;
		ProgramGlobals::DirectionEnum direction;
		bool                          debug;
		bool                          enablePersistentSvd;
		bool                          serialSvd;
	};

	using BuildingBlockType = typename BlockDiagonalMatrixType::BuildingBlockType;

	virtual ~DensityMatrixBase() { }

	virtual const BlockDiagonalMatrixType& operator()() = 0;

	virtual void diag(typename PsimagLite::Vector<RealType>::Type&, char) = 0;

	virtual const typename PsimagLite::Vector<MatrixType>::Type& vts() const
	{
		return vtsEmpty_;
	}

	virtual const typename PsimagLite::Vector<VectorRealType>::Type& s() const
	{
		return sEmpty_;
	}

	virtual const typename BasisWithOperatorsType::VectorQnType& qns() const
	{
		return qnsEmpty_;
	}

private:

	typename PsimagLite::Vector<MatrixType>::Type     vtsEmpty_;
	typename PsimagLite::Vector<VectorRealType>::Type sEmpty_;
	typename BasisWithOperatorsType::VectorQnType     qnsEmpty_;

}; // class DensityMatrixBase
} // namespace Dmrg

/*@}*/
#endif
