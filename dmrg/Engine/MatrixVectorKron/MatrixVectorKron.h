// SPDX-FileCopyrightText: Copyright (c) 2009, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 2.0.0]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file MatrixVectorKron.h
 *
 *  A class to encapsulate the product x+=Hy,
 *  where x and y are vectors and H is the Hamiltonian matrix
 *
 */
#ifndef MATRIX_VECTOR_KRON_H
#define MATRIX_VECTOR_KRON_H

#include "InitKronHamiltonian.h"
#include "KronMatrix.h"
#include "MatrixVectorBase.h"
#include <PsimagLite/Vector.h>

namespace Dmrg {
template <typename ComplexOrRealType_>
class MatrixVectorKron final : public MatrixVectorBase<ComplexOrRealType_> {

	using BaseType = MatrixVectorBase<ComplexOrRealType_>;

	static const bool CHECK_KRON = true;

public:

	using ModelType                 = typename BaseType::ModelType;
	using ModelHelperType           = typename ModelType::ModelHelperType;
	using ParametersType            = typename ModelType::ParametersType;
	using RealType                  = typename ModelHelperType::RealType;
	using InitKronType              = InitKronHamiltonian<ModelType>;
	using KronMatrixType            = KronMatrix<InitKronType>;
	using SparseMatrixType          = typename ModelHelperType::SparseMatrixType;
	using ComplexOrRealType         = ComplexOrRealType_;
	using VectorRealType            = typename PsimagLite::Vector<RealType>::Type;
	using VectorType                = typename PsimagLite::Vector<ComplexOrRealType>::Type;
	using FullMatrixType            = PsimagLite::Matrix<ComplexOrRealType>;
	using value_type                = typename SparseMatrixType::value_type;
	using HamiltonianConnectionType = typename ModelType::HamiltonianConnectionType;
	using MatrixSolverEnum          = typename ParametersType::MatrixSolverEnum;

	MatrixVectorKron(const ModelType&                     model,
	                 const HamiltonianConnectionType&     hc,
	                 const typename ModelHelperType::Aux& aux)
	    : BaseType(hc, aux)
	    , params_(model.params())
	    , isLdaggerL_(validatedLdaggerL(model))
	    , initKron_(model, hc, aux)
	    , kronMatrix_(initKron_, "Hamiltonian")
	    , time_(0, 0)
	{
		int maxMatrixRankStored = model.params().maxMatrixRankStored;
		if (hc.modelHelper().size(aux.m()) > maxMatrixRankStored)
			return;

		hc.fullHamiltonian(matrixStored_, aux, model.isHermitian());
		if (isLdaggerL_)
			transposeConjugate(transpose_, matrixStored_);
		else
			assert(isHermitian(matrixStored_, true));

		if (!isLdaggerL_)
			checkKron();
	}

	~MatrixVectorKron() override
	{
		std::cout << "DeltaClock matrixVectorProduct " << time_.millis() << "\n";
	}

	void matrixVectorProduct(VectorType& x, const VectorType& y) const override
	{
		const PsimagLite::MemoryUsage::TimeHandle time1
		    = PsimagLite::ProgressIndicator::time();

		if (!isLdaggerL_) {
			if (matrixStored_.rows() > 0)
				matrixStored_.matrixVectorProduct(x, y);
			else
				kronMatrix_.matrixVectorProduct(x, y);
		} else if (matrixStored_.rows() > 0) {
			VectorType intermediate(y.size(), 0.0);
			matrixStored_.matrixVectorProduct(intermediate, y);
			transpose_.matrixVectorProduct(x, intermediate);
		} else {
			VectorType intermediate(y.size(), 0.0);
			kronMatrix_.matrixVectorProduct(intermediate, y);
			kronMatrix_.matrixVectorProductAdjoint(x, intermediate);
		}

		const PsimagLite::MemoryUsage::TimeHandle time2
		    = PsimagLite::ProgressIndicator::time();
		const PsimagLite::MemoryUsage::TimeHandle deltaTime = time2 - time1;
		time_ += deltaTime;
	}

	void fullDiag(VectorRealType& eigs, FullMatrixType& fm) const override
	{
		if (isLdaggerL_)
			throw PsimagLite::RuntimeError(
			    "MatrixVectorKron::fullDiag is unavailable with LdaggerL\n");

		BaseType::fullDiag(eigs, fm, matrixStored_, params_.maxMatrixRankStored);
	}

private:

	static bool validatedLdaggerL(const ModelType& model)
	{
		if (!model.params().options.isSet("LdaggerL"))
			return false;

		switch (model.params().matrix_solver_enum) {
		case MatrixSolverEnum::LANCZOS:
			return true;
		case MatrixSolverEnum::DENSE:
			throw PsimagLite::RuntimeError(
			    "LdaggerL with MatrixVectorKron does not support MatrixSolver=Dense\n");
		case MatrixSolverEnum::ARNOLDISAI:
			throw PsimagLite::RuntimeError(
			    "LdaggerL with MatrixVectorKron does not support MatrixSolver=ArnoldiSaI\n");
		}

		throw PsimagLite::RuntimeError(
		    "LdaggerL with MatrixVectorKron supports only MatrixSolver=Lanczos\n");
	}

	void checkKron() const
	{
		if (!CHECK_KRON)
			return;

#ifdef NDEBUG
		return;
#endif

		SizeType n = this->rows();
		std::cout << n << "\n";
		FullMatrixType m(n, n);
		for (SizeType i = 0; i < n; ++i) {
			VectorType e(n, 0.0);
			e[i] = 1.0;
			VectorType ey(n, 0.0);
			kronMatrix_.matrixVectorProduct(ey, e);
			for (SizeType j = 0; j < n; ++j)
				m(i, j) = ey[j];
		}

		std::cout << "Matrix as thought of by Kron\n";
		std::cout << m;
		assert(isHermitian(m));
		std::cout << "Correct matrix\n";
		std::cout << matrixStored_;
	}

	const ParametersType&                       params_;
	const bool                                  isLdaggerL_;
	InitKronType                                initKron_;
	KronMatrixType                              kronMatrix_;
	SparseMatrixType                            matrixStored_;
	SparseMatrixType                            transpose_;
	mutable PsimagLite::MemoryUsage::TimeHandle time_;
}; // class MatrixVectorKron
} // namespace Dmrg

/*@}*/
#endif
