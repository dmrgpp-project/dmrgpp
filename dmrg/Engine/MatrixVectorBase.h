// SPDX-FileCopyrightText: Copyright (c) 2009-2014, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file MatrixVectorBase.h
 *
 *  The base of a hierarchy to encapsulate the product x+=Hy,
 *  where x and y are vectors and H is the Hamiltonian matrix
 *
 */
#ifndef DMRG_MATRIX_VECTOR_BASE_H
#define DMRG_MATRIX_VECTOR_BASE_H

#include "MatrixVectorTypes.hpp"
#include <PsimagLite/PsimagLite.h>
#include <cassert>
#include <vector>

namespace Dmrg {
template <typename ComplexOrRealType_, typename TypesType_ = MatrixVectorTypes<ComplexOrRealType_>>
class MatrixVectorBase {

public:

	using TypesType         = TypesType_;
	using ModelType         = typename TypesType::ModelType;
	using ModelHelperType   = typename TypesType::ModelHelperType;
	using RealType          = typename TypesType::RealType;
	using SparseMatrixType  = typename TypesType::SparseMatrixType;
	using ComplexOrRealType = ComplexOrRealType_;
	using value_type        = ComplexOrRealType;
	using VectorRealType    = typename PsimagLite::Vector<RealType>::Type;
	using VectorType        = typename PsimagLite::Vector<ComplexOrRealType>::Type;
	using FullMatrixType    = PsimagLite::Matrix<ComplexOrRealType>;

	using HamiltonianConnectionType = typename ModelType::HamiltonianConnectionType;
	using AuxType                   = typename ModelHelperType::Aux;

	MatrixVectorBase(const HamiltonianConnectionType& hc, const AuxType& aux)
	    : rows_(computeRows(hc, aux))
	{ }

	virtual ~MatrixVectorBase() = default;

	SizeType rows() const { return rows_; }

	SizeType cols() const { return rows_; }

	virtual void matrixVectorProduct(VectorType& x, const VectorType& y) const = 0;

	SizeType reflectionSector() const { return 0; }

	void reflectionSector(SizeType) { }

	virtual void fullDiag(VectorRealType& eigs, FullMatrixType& fm) const = 0;

	virtual const SparseMatrixType& toCRS() const
	{
		throw PsimagLite::RuntimeError("This MatrixVector class doesn't support toCRS()\n");
	}

	static void fullDiag(VectorRealType&         eigs,
	                     FullMatrixType&         fm,
	                     const SparseMatrixType& matrixStored,
	                     int                     tmp)
	{
		SizeType maxMatrixRankStored = (tmp < 0) ? 0 : tmp;

		if (matrixStored.rows() == 0) {
			PsimagLite::String str("MatrixVectorBase:fullDiag: no stored matrix\n");
			str += "\trow= " + ttos(eigs.size()) + " max row= ";
			str += ttos(maxMatrixRankStored) + "\n";
			str += "Please add or increase MaxMatrixRankStored="
			    + ttos(2 + matrixStored.rows());
			str += " in your input file\n";
			err(str);
		}

		if (matrixStored.rows() > maxMatrixRankStored) {
			PsimagLite::String str("MatrixVectorBase:fullDiag: internal error!\n");
			str += "\trow= " + ttos(matrixStored.rows()) + " max row= ";
			str += ttos(maxMatrixRankStored) + "\n";
			str += "Please add or increase MaxMatrixRankStored= to at least "
			    + ttos(2 + eigs.size()) + " in your input file\n";
			err(str);
		}

		fm = matrixStored.toDense();
		diag(fm, eigs, 'V');
	}

private:

	static SizeType computeRows(const HamiltonianConnectionType& hc, const AuxType& aux)
	{
		const auto&    lrs    = hc.modelHelper().leftRightSuper();
		const SizeType m      = aux.m();
		const SizeType offset = lrs.super().partition(m);
		assert(lrs.super().partition(m + 1) >= offset);
		const SizeType rows = lrs.super().partition(m + 1) - offset;
		assert(rows > 0);
		return rows;
	}

	const SizeType rows_;
}; // class MatrixVectorBase
} // namespace Dmrg

/*@}*/
#endif
