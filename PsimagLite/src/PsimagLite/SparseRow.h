// SPDX-FileCopyrightText: Copyright (c) 2009, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file SparseRow.h
 *
 *  Sparse row allows one to accumulate non-zeros to a row
 *  of a CrsMatrix and add the row at the end
 */
#ifndef SPARSE_ROW_H
#define SPARSE_ROW_H

#include "Sort.h"
#include <cassert>

namespace PsimagLite {

template <typename CrsMatrixType> class SparseRow {
public:

	using ValueType   = typename CrsMatrixType::value_type;
	using ColumnsType = typename Vector<SizeType>::Type;
	using VectorType  = typename Vector<ValueType>::Type;

	void add(SizeType col, ValueType value)
	{
		cols_.push_back(col);
		values_.push_back(value);
	}

	ValueType matrixVectorProduct(const VectorType& y) const
	{
		ValueType sum = 0;
		for (SizeType i = 0; i < cols_.size(); i++)
			sum += values_[i] * y[cols_[i]];
		return sum;
	}
	//		void clear()
	//		{
	//			cols_.clear();
	//			values_.clear();
	//		}

	SizeType finalize(CrsMatrixType& matrix)
	{
		assert(cols_.size() == values_.size());
		if (cols_.size() == 0)
			return 0;

		Sort<ColumnsType> s;
		ColumnsType       iperm(cols_.size());
		s.sort(cols_, iperm);
		SizeType  prevCol = cols_[0];
		SizeType  counter = 0;
		ValueType value   = 0;
		for (SizeType i = 0; i < cols_.size(); i++) {
			assert(cols_[i] < matrix.cols());
			if (cols_[i] == prevCol) {
				value += values_[iperm[i]];
				continue;
			}
			matrix.pushCol(prevCol);
			matrix.pushValue(value);
			counter++;
			value   = values_[iperm[i]];
			prevCol = cols_[i];
		}
		matrix.pushCol(prevCol);
		matrix.pushValue(value);
		counter++;
		return counter;
	}

	ValueType finalize(const VectorType& y)
	{
		assert(cols_.size() == values_.size());
		if (cols_.size() == 0)
			return 0;

		ValueType sum = 0.0;
		for (SizeType i = 0; i < cols_.size(); i++)
			sum += values_[i] * y[cols_[i]];
		return sum;
	}

private:

	ColumnsType cols_;
	VectorType  values_;

}; // class SparseRow

} // namespace PsimagLite

/*@}*/
#endif // SPARSE_ROW_H
