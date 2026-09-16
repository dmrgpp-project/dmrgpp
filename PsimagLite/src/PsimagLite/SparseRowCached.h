// SPDX-FileCopyrightText: Copyright (c) 2009, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file SparseRowCached.h
 *
 *  Sparse row allows one to accumulate non-zeros to a row
 *  of a CrsMatrix and add the row at the end
 */
#ifndef SPARSE_ROW_CACHED_H
#define SPARSE_ROW_CACHED_H

#include "Sort.h"
#include <cassert>

namespace PsimagLite {

template <typename CrsMatrixType> class SparseRowCached {

public:

	using ValueType   = typename CrsMatrixType::value_type;
	using ColumnsType = typename Vector<SizeType>::Type;
	using VectorType  = typename Vector<ValueType>::Type;

	SparseRowCached(SizeType cacheSize)
	    : cols_(cacheSize)
	    , values_(cacheSize)
	    , counter_(0)
	{ }

	void add(SizeType col, ValueType value)
	{
		/*			cols_.push_back(col);
		                        values_.push_back(value);*/
		SizeType cacheSize = cols_.size();
		if (counter_ >= cacheSize) {
			cacheSize *= 2;
			cols_.resize(cacheSize);
			values_.resize(cacheSize);
		}

		cols_[counter_]   = col;
		values_[counter_] = value;
		counter_++;
	}

	ValueType matrixVectorProduct(const VectorType& y)
	{
		ValueType sum = 0;
		for (SizeType i = 0; i < counter_; i++)
			sum += values_[i] * y[cols_[i]];
		counter_ = 0;
		return sum;
	}
	//		void clear()
	//		{
	//			cols_.clear();
	//			values_.clear();
	//		}

	SizeType finalize(CrsMatrixType& matrix)
	{
		// assert(cols_.size()==values_.size());
		if (counter_ == 0)
			return 0;

		Sort<ColumnsType> s;
		ColumnsType       iperm(counter_);
		s.sort(cols_, iperm, counter_);
		SizeType  prevCol = cols_[0];
		SizeType  counter = 0;
		ValueType value   = 0;
		for (SizeType i = 0; i < counter_; i++) {
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
		counter_ = 0;
		return counter;
	}

private:

	ColumnsType                      cols_;
	typename Vector<ValueType>::Type values_;
	SizeType                         counter_;

}; // class SparseRowCached

} // namespace PsimagLite

/*@}*/
#endif // SPARSE_ROW_CACHED_H
