// SPDX-FileCopyrightText: Copyright (c) 2009-2016, 2017, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [Lanczos, Version 2.]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup LanczosPlusPlus */
/*@{*/

/*! \file InternalProductOnTheFly.h
 *
 *  A class to encapsulate the product x+=Hy, where x and y are vectors and H is the Hamiltonian
 * matrix
 *
 */
#ifndef INTERNALPRODUCT_OTF_H
#define INTERNALPRODUCT_OTF_H

#include <PsimagLite/Matrix.h>
#include <PsimagLite/Vector.h>
#include <cassert>
#include <vector>

namespace LanczosPlusPlus {
template <typename ModelType_, typename SpecialSymmetryType_> class InternalProductOnTheFly {

public:

	typedef ModelType_                                           ModelType;
	typedef SpecialSymmetryType_                                 SpecialSymmetryType;
	typedef typename ModelType::BasisBaseType                    BasisType;
	typedef typename SpecialSymmetryType::SparseMatrixType       SparseMatrixType;
	typedef typename ModelType::RealType                         RealType;
	typedef typename ModelType::GeometryType                     GeometryType;
	typedef typename GeometryType::ComplexOrRealType             ComplexOrRealType;
	typedef ComplexOrRealType                                    value_type;
	typedef PsimagLite::Matrix<ComplexOrRealType>                MatrixType;
	typedef typename PsimagLite::Vector<RealType>::Type          VectorRealType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;

	InternalProductOnTheFly(const ModelType& model,
	                        const BasisType& basis,
	                        SpecialSymmetryType&)
	    : model_(model)
	    , basis_(basis)
	{ }

	InternalProductOnTheFly(const ModelType& model, SpecialSymmetryType&)
	    : model_(model)
	    , basis_(model.basis())
	{ }

	SizeType rows() const { return basis_.size(); }

	void matrixVectorProduct(VectorType& x, const VectorType& y) const
	{
		model_.matrixVectorProduct(x, y, basis_);
	}

	SizeType reflectionSector() const { return 0; }

	void specialSymmetrySector(SizeType) { }

	void fullDiag(VectorRealType&, MatrixType&)
	{
		err("no fullDiag possible when on the fly\n");
	}

private:

	const ModelType& model_;
	const BasisType& basis_;
}; // class InternalProductOnTheFly
} // namespace LanczosPlusPlus

/*@}*/
#endif
