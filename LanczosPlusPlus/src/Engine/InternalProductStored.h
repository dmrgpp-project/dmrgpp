// SPDX-FileCopyrightText: Copyright (c) 2010-2014, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [Lanczos++, Version 1.0.0]
// [by G.A., Oak Ridge National Laboratory]


/** \ingroup LanczosPlusPlus */
/*@{*/

/*! \file InternalProductStored.h
 *
 *  A class to encapsulate the product x+=Hy,
 *  where x and y are vectors and H is the Hamiltonian matrix
 *
 */
#ifndef InternalProductStored_HEADER_H
#define InternalProductStored_HEADER_H

#include <PsimagLite/Vector.h>
#include <vector>

namespace LanczosPlusPlus {
template <typename ModelType_, typename SpecialSymmetryType_> class InternalProductStored {

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

	InternalProductStored(const ModelType&     model,
	                      const BasisType&     basis,
	                      SpecialSymmetryType& rs)
	    : rs_(rs)
	{
		rs_.init(model, basis);
	}

	InternalProductStored(const ModelType& model, SpecialSymmetryType& rs)
	    : rs_(rs)
	{
		rs_.init(model, model.basis());
	}

	SizeType rows() const { return rs_.rows(); }

	void matrixVectorProduct(VectorType& x, const VectorType& y) const
	{
		rs_.matrixVectorProduct(x, y);
	}

	void specialSymmetrySector(SizeType p) { rs_.setPointer(p); }

	void fullDiag(VectorRealType& eigs, MatrixType& z) { rs_.fullDiag(eigs, z); }

private:

	SpecialSymmetryType& rs_;
}; // class InternalProductStored
} // namespace LanczosPlusPlus

/*@}*/
#endif
