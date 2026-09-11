/*
Copyright (c) 2009, 2026, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 6+]
[by G.A., Oak Ridge National Laboratory]

UT Battelle Open Source Software License 11242008

OPEN SOURCE LICENSE

Subject to the conditions of this License, each
contributor to this software hereby grants, free of
charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), a
perpetual, worldwide, non-exclusive, no-charge,
royalty-free, irrevocable copyright license to use, copy,
modify, merge, publish, distribute, and/or sublicense
copies of the Software.

1. Redistributions of Software must retain the above
copyright and license notices, this list of conditions,
and the following disclaimer.  Changes or modifications
to, or derivative works of, the Software should be noted
with comments and the contributor and organization's
name.

2. Neither the names of UT-Battelle, LLC or the
Department of Energy nor the names of the Software
contributors may be used to endorse or promote products
derived from this software without specific prior written
permission of UT-Battelle.

3. The software and the end-user documentation included
with the redistribution, with or without modification,
must include the following acknowledgment:

"This product includes software produced by UT-Battelle,
LLC under Contract No. DE-AC05-00OR22725  with the
Department of Energy."

*********************************************************
DISCLAIMER

THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT OWNER, CONTRIBUTORS, UNITED STATES GOVERNMENT,
OR THE UNITED STATES DEPARTMENT OF ENERGY BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED
STATES DEPARTMENT OF ENERGY, NOR THE COPYRIGHT OWNER, NOR
ANY OF THEIR EMPLOYEES, REPRESENTS THAT THE USE OF ANY
INFORMATION, DATA, APPARATUS, PRODUCT, OR PROCESS
DISCLOSED WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*********************************************************

*/
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
