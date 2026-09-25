// SPDX-FileCopyrightText: Copyright (c) 2009-2015, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file MatrixOrIdentity.h
 *
 *
 *
 */

#ifndef MATRIX_OR_IDENTITY
#define MATRIX_OR_IDENTITY

namespace Dmrg {

template <typename SparseMatrixType> class MatrixOrIdentity {

	using SparseElementType = typename SparseMatrixType::value_type;

public:

	MatrixOrIdentity(bool enabled, const SparseMatrixType& m)
	    : enabled_(enabled)
	    , m_(m)
	    , one_(1.0)
	{ }

	SizeType getRowPtr(SizeType i) const { return (enabled_) ? m_.getRowPtr(i) : i; }

	int getColOrExit(SizeType i) const
	{
		if (enabled_)
			return m_.getCol(i);
		if (i < m_.cols())
			return i;
		return -1;
	}

	const SparseElementType& getValue(SizeType i) const
	{
		return (enabled_) ? m_.getValue(i) : one_;
	}

private:

	bool                    enabled_;
	const SparseMatrixType& m_;
	SparseElementType       one_;
}; // class MatrixOrIdentity

} // namespace Dmrg
/*@}*/
#endif // MATRIX_OR_IDENTITY
