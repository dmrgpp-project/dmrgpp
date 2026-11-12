
// SPDX-FileCopyrightText: Copyright (c) 2009, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 0.0.1]
// [DMRG++ authors; see AUTHORS.]

#ifndef SAMPLE_CRSMATRIX_HEADER_H
#define SAMPLE_CRSMATRIX_HEADER_H

#include "Sort.h"
#include "loki/TypeTraits.h"
#include <cassert>
#include <cstdlib> // <-- just for drand48 and srand48(seed)
#include <iostream>
#include <stdexcept>
#include <vector>

namespace PsimagLite {

template <typename T> class SampleCRSMatrix {

public:

	using value_type = T;

	SampleCRSMatrix(SizeType rank)
	    : rank_(rank)
	    , rowptr_(rank + 1)
	{ }

	SampleCRSMatrix(SizeType rank, SizeType seed, SizeType nonZeros, T maxValue)
	    : rank_(rank)
	    , rowptr_(rank + 1)
	{
		srand48(seed);
		typename Vector<SizeType>::Type rows, cols;
		typename Vector<T>::Type        vals;
		for (SizeType i = 0; i < nonZeros; i++) {
			// pick a row
			SizeType row = SizeType(drand48() * rank);
			// and a column
			SizeType col = SizeType(drand48() * rank);
			// and a value
			T val = drand48() * maxValue;
			rows.push_back(row);
			cols.push_back(col);
			vals.push_back(val);
		}
		// fill the matrix with this data:
		fillMatrix(rows, cols, vals);
	}

	template <typename SomeIoInputType> SampleCRSMatrix(SomeIoInputType& io)
	{
		io >> rank_;
		readVector(io, rowptr_);
		readVector(io, colind_);
		readVector(io, values_);
	}

	void setRow(int n, int v) { rowptr_[n] = v; }

	void pushCol(int i) { colind_.push_back(i); }

	void pushValue(const T& value) { values_.push_back(value); }

	void matrixVectorProduct(typename Vector<T>::Type&       x,
	                         const typename Vector<T>::Type& y) const
	{
		for (SizeType i = 0; i < y.size(); i++)
			for (SizeType j = rowptr_[i]; j < rowptr_[i + 1]; j++)
				x[i] += values_[j] * y[colind_[j]];
	}

	SizeType rank() const { return rank_; }

	template <typename SomeIoOutputType> void save(SomeIoOutputType& io) const
	{
		io << rank_ << "\n";
		saveVector(io, rowptr_);
		saveVector(io, colind_);
		saveVector(io, values_);
	}

	SizeType getRowPtr(SizeType i) const
	{
		assert(i < rowptr_.size());
		return rowptr_[i];
	}

	SizeType getCol(SizeType i) const
	{
		assert(i < colind_.size());
		return colind_[i];
	}

	const T& getValues(SizeType i) const
	{
		assert(i < values_.size());
		return values_[i];
	}

private:

	template <typename SomeIoOutputType, typename SomeVectorType>
	typename EnableIf<IsVectorLike<SomeVectorType>::True, void>::Type
	saveVector(SomeIoOutputType& io, const SomeVectorType& v) const
	{
		io << v.size() << "\n";
		for (SizeType i = 0; i < v.size(); i++) {
			io << v[i] << " ";
		}
		io << "\n";
	}

	template <typename SomeIoInputType, typename SomeVectorType>
	typename EnableIf<IsVectorLike<SomeVectorType>::True, void>::Type
	readVector(SomeIoInputType& io, SomeVectorType& v) const
	{
		int size = 0;
		io >> size;
		if (size < 0)
			throw RuntimeError("readVector: size is zero\n");
		v.resize(size);
		for (SizeType i = 0; i < v.size(); i++) {
			io >> v[i];
		}
	}

	void fillMatrix(typename Vector<SizeType>::Type& rows,
	                typename Vector<SizeType>::Type& cols,
	                typename Vector<T>::Type&        vals)
	{
		Sort<typename Vector<SizeType>::Type> s;
		typename Vector<SizeType>::Type       iperm(rows.size());
		s.sort(rows, iperm);
		SizeType counter = 0;
		SizeType prevRow = rows[0] + 1;
		for (SizeType i = 0; i < rows.size(); i++) {
			SizeType row = rows[i];
			if (prevRow != row) {
				// add new row
				rowptr_[row] = counter++;
				prevRow      = row;
			}
			colind_.push_back(cols[iperm[i]]);
			values_.push_back(vals[iperm[i]]);
		}
		SizeType lastNonZeroRow = rows[rows.size() - 1];
		for (SizeType i = lastNonZeroRow + 1; i <= rank_; i++)
			rowptr_[i] = counter;
	}

	SizeType                        rank_;
	typename Vector<SizeType>::Type rowptr_;
	typename Vector<SizeType>::Type colind_;
	typename Vector<T>::Type        values_;

}; // class SampleCRSMatrix

template <typename T> std::ostream& operator<<(std::ostream& os, const SampleCRSMatrix<T>& m)
{
	m.save(os);
	return os;
}

} // namespace PsimagLite
#endif
