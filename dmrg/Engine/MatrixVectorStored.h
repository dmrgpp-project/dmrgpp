// SPDX-FileCopyrightText: Copyright (c) 2009, 2026, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 6+]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file MatrixVectorStored.h
 *
 *  A class to encapsulate the product x+=Hy,
 *  where x and y are vectors and H is the Hamiltonian matrix
 *
 */
#ifndef MatrixVectorStored_HEADER_H
#define MatrixVectorStored_HEADER_H

#include "MatrixVectorBase.h"
#include <PsimagLite/ProgressIndicator.h>
#include <vector>

namespace Dmrg {
template <typename ComplexOrRealType_, typename TypesType_ = MatrixVectorTypes<ComplexOrRealType_>>
class MatrixVectorStored final : public MatrixVectorBase<ComplexOrRealType_, TypesType_> {

	using BaseType = MatrixVectorBase<ComplexOrRealType_, TypesType_>;

public:

	using ModelType                 = typename BaseType::ModelType;
	using HamiltonianConnectionType = typename ModelType::HamiltonianConnectionType;
	using ParametersType            = typename ModelType::ParametersType;
	using ModelHelperType           = typename ModelType::ModelHelperType;
	using SparseMatrixType          = typename ModelHelperType::SparseMatrixType;
	using RealType                  = typename ModelHelperType::RealType;
	using value_type                = typename SparseMatrixType::value_type;
	using ComplexOrRealType         = ComplexOrRealType_;
	using VectorRealType            = typename PsimagLite::Vector<RealType>::Type;
	using VectorType                = typename BaseType::VectorType;
	using OptionsType               = typename ParametersType::OptionsType;
	using MatrixSolverEnum          = typename ParametersType::MatrixSolverEnum;
	using FullMatrixType            = PsimagLite::Matrix<ComplexOrRealType>;

	MatrixVectorStored(const ModelType&                     model,
	                   const HamiltonianConnectionType&     hc,
	                   const typename ModelHelperType::Aux& aux)
	    : BaseType(hc, aux)
	    , model_(model)
	    , isLdaggerL_(validatedLdaggerL(model))
	    , progress_("MatrixVectorStored")
	{
		const OptionsType& options     = model.params().options;
		const bool         debugMatrix = options.isSet("debugmatrix");

		matrixStored_.clear();

		hc.fullHamiltonian(matrixStored_, aux, model.isHermitian());
		if (isLdaggerL_) {
			transposeConjugate(transpose_, matrixStored_);
		}

		PsimagLite::OstringStream                     msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg  = msgg();
		SizeType                                      rows = matrixStored_.rows();
		msg << "fullHamiltonian has rank=" << rows;
		msg << " nonzeros=" << matrixStored_.nonZeros();
		progress_.printline(msgg, std::cout);
		if (debugMatrix)
			printFullMatrix(matrixStored_, "matrix", 1);
		if (rows > 5000)
			std::cerr << "WARNING: MatrixVectorStored being used for a large run!\n";
	}

	const SparseMatrixType& toCRS() const override
	{
		if (isLdaggerL_)
			throw PsimagLite::RuntimeError(
			    "MatrixVectorStored::toCRS is unavailable with LdaggerL\n");

		return matrixStored_;
	}

	void matrixVectorProduct(VectorType& x, const VectorType& y) const override
	{
		if (!isLdaggerL_) {
			matrixStored_.matrixVectorProduct(x, y);
			return;
		}

		VectorType intermediate(y.size(), 0.0);
		matrixStored_.matrixVectorProduct(intermediate, y);
		transpose_.matrixVectorProduct(x, intermediate);
	}

	void fullDiag(VectorRealType& eigs, FullMatrixType& fm) const override
	{
		if (isLdaggerL_)
			throw PsimagLite::RuntimeError(
			    "MatrixVectorStored::fullDiag is unavailable with LdaggerL\n");

		BaseType::fullDiag(eigs, fm, matrixStored_, model_.params().maxMatrixRankStored);
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
			    "LdaggerL with MatrixVectorStored does not support MatrixSolver=Dense\n");
		case MatrixSolverEnum::ARNOLDISAI:
			throw PsimagLite::RuntimeError(
			    "LdaggerL with MatrixVectorStored does not support MatrixSolver=ArnoldiSaI\n");
		}

		throw PsimagLite::RuntimeError(
		    "LdaggerL with MatrixVectorStored supports only MatrixSolver=Lanczos\n");
	}

	const ModelType&              model_;
	const bool                    isLdaggerL_;
	SparseMatrixType              matrixStored_;
	SparseMatrixType              transpose_;
	PsimagLite::ProgressIndicator progress_;
}; // class MatrixVectorStored
} // namespace Dmrg

/*@}*/
#endif
