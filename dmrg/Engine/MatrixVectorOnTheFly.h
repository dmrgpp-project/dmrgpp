// SPDX-FileCopyrightText: Copyright (c) 2009, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 2.0.0]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file MatrixVectorOnTheFly.h
 *
 *  A class to encapsulate the product x+=Hy,
 *  where x and y are vectors and H is the Hamiltonian matrix
 *
 */
#ifndef MATRIX_VECTOR_OTF_H
#define MATRIX_VECTOR_OTF_H

#include "MatrixVectorBase.h"
#include <vector>

namespace Dmrg {
template <typename ComplexOrRealType_>
class MatrixVectorOnTheFly final : public MatrixVectorBase<ComplexOrRealType_> {

	using BaseType = MatrixVectorBase<ComplexOrRealType_>;

public:

	using ModelType                 = typename BaseType::ModelType;
	using ModelHelperType           = typename ModelType::ModelHelperType;
	using RealType                  = typename ModelHelperType::RealType;
	using SparseMatrixType          = typename ModelHelperType::SparseMatrixType;
	using value_type                = typename SparseMatrixType::value_type;
	using ComplexOrRealType         = ComplexOrRealType_;
	using VectorRealType            = typename PsimagLite::Vector<RealType>::Type;
	using VectorType                = typename BaseType::VectorType;
	using FullMatrixType            = PsimagLite::Matrix<ComplexOrRealType>;
	using HamiltonianConnectionType = typename ModelType::HamiltonianConnectionType;
	using AuxType                   = typename ModelHelperType::Aux;

	MatrixVectorOnTheFly(const ModelType&                 model,
	                     const HamiltonianConnectionType& hc,
	                     const AuxType&                   aux)
	    : BaseType(hc, aux)
	    , model_(model)
	    , hc_(hc)
	    , aux_(aux)
	{
		if (model.params().options.isSet("LdaggerL"))
			throw PsimagLite::RuntimeError(
			    "LdaggerL is not implemented for MatrixVectorOnTheFly\n");

		int maxMatrixRankStored = model.params().maxMatrixRankStored;
		if (hc.modelHelper().size(aux_.m()) > maxMatrixRankStored)
			return;

		hc.fullHamiltonian(matrixStored_, aux_, model_.isHermitian());
		assert(isHermitian(matrixStored_, true));
	}

	void matrixVectorProduct(VectorType& x, const VectorType& y) const override
	{
		if (matrixStored_.rows() > 0)
			matrixStored_.matrixVectorProduct(x, y);
		else
			model_.matrixVectorProduct(x, y, hc_, aux_);
	}

	void fullDiag(VectorRealType& eigs, FullMatrixType& fm) const override
	{
		int mrs = model_.params().maxMatrixRankStored;
		if (mrs < static_cast<int>(this->rows())) {
			std::cerr << "Full diag will likely fail, it would need ";
			std::cerr << this->rows() << " but you gave only " << mrs << "\n";
		}

		BaseType::fullDiag(eigs, fm, matrixStored_, mrs);
	}

private:

	const ModelType&                 model_;
	const HamiltonianConnectionType& hc_;
	const AuxType&                   aux_;
	SparseMatrixType                 matrixStored_;
}; // class MatrixVectorOnTheFly
} // namespace Dmrg

/*@}*/
#endif
