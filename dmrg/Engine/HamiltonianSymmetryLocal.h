// SPDX-FileCopyrightText: Copyright (c) 2009, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 2.0.0]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file HamiltonianSymmetryLocal.h
 *
 *  This class contains the implementation of local symmetries
 *  An object of this class is contained by DmrgBasisImplementation
 *
 */
#ifndef HAM_SYMM_LOCAL_H
#define HAM_SYMM_LOCAL_H
#include "ProgramGlobals.h"
#include <PsimagLite/Io/IoSelector.h>
#include <PsimagLite/Sort.h>

namespace Dmrg {
template <typename SparseMatrixType> class HamiltonianSymmetryLocal {

	using SparseElementType = typename SparseMatrixType::value_type;
	using RealType          = typename PsimagLite::Real<SparseElementType>::Type;
	using FactorsType       = PsimagLite::CrsMatrix<RealType>;
	using VectorSizeType    = typename PsimagLite::Vector<SizeType>::Type;

public:

	SizeType getFlavor(SizeType) const
	{
		return 0; // meaningless
	}

	template <typename SolverParametersType>
	void calcRemovedIndices(VectorSizeType&       removedIndices,
	                        const VectorSizeType& perm,
	                        SizeType              kept,
	                        const SolverParametersType&) const
	{
		const SizeType permSize = perm.size();
		if (permSize <= kept)
			return;

		SizeType target = permSize - kept;

		removedIndices.clear();
		for (SizeType i = 0; i < target; ++i) {
			if (removedIndices.size() >= target)
				break;
			if (PsimagLite::indexOrMinusOne(removedIndices, perm[i]) >= 0)
				continue;
			removedIndices.push_back(perm[i]);
		}
	}

	template <typename IoInputter>
	void
	read(IoInputter&,
	     PsimagLite::String,
	     bool,
	     typename PsimagLite::EnableIf<PsimagLite::IsInputLike<IoInputter>::True, int>::Type
	     = 0)
	{ }

	void write(PsimagLite::IoSelector::Out&          io,
	           PsimagLite::String                    label,
	           PsimagLite::IoNgSerializer::WriteMode mode) const
	{
		io.write(0, label + "FACTORSSIZE", mode);
	}
}; // class HamiltonianSymmetryLocal
} // namespace Dmrg

/*@}*/
#endif
