// SPDX-FileCopyrightText: Copyright (c) 2009-2013, 2026, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup PsimagLite */
/*@{*/

/*! \file MatrixSolverBase.h
 *
 *  Virtual class to select LanczosSolver or ArnoldiSaI
 *
 */

#ifndef PSI_MATRIX_SOLVER_BASE_H
#define PSI_MATRIX_SOLVER_BASE_H
#include "Complex.h"
#include <cassert>
#include <vector>

namespace PsimagLite {

template <typename MatrixType_> class MatrixSolverBase {

public:

	using MatrixType        = MatrixType_;
	using ComplexOrRealType = typename MatrixType_::value_type;
	using VectorType        = std::vector<ComplexOrRealType>;
	using RealType          = typename Real<ComplexOrRealType>::Type;
	using VectorRealType    = std::vector<RealType>;
	using VectorVectorType  = std::vector<VectorType>;

	// To avoid compiler warnings
	virtual ~MatrixSolverBase() = default;

	virtual void computeOneState(RealType&, VectorType&, const VectorType&, SizeType) = 0;

	virtual void
	computeAllStatesBelow(VectorRealType&, VectorVectorType&, const VectorType&, SizeType)
	    = 0;

	static bool isReorthoEnabled(const std::string& options, bool lotaMemory)
	{
		if (options.find("reortho") == PsimagLite::String::npos)
			return false;

		bool canReortho = (lotaMemory);

		if (!canReortho) {
			throw RuntimeError(
			    "MatrixSolverBase: Reortho requested but cannot. Suggestion: Delete "
			    "reortho from input or set lotaMemory=true\n");
		}

		return true;
	}

}; // class MatrixSolverBase

} // namespace PsimagLite

/*@}*/
#endif // PSI_MATRIX_SOLVER_BASE_H
