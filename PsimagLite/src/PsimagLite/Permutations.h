// SPDX-FileCopyrightText: Copyright (c) 2009-2012, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file Permutations.h
 *
 * Raw computations for a free Hubbard model
 *
 */
#ifndef PERMUTATIONS_H_H
#define PERMUTATIONS_H_H

#include "Complex.h" // in PsimagLite
#include "Sort.h"

namespace PsimagLite {

template <typename ContainerType> class Permutations {
	using FieldType = typename ContainerType::value_type;

public:

	using value_type = FieldType;

	Permutations(const ContainerType& orig)
	    : data_(orig.size())
	{
		for (SizeType i = 0; i < data_.size(); i++)
			data_[i] = orig[i];
		Sort<typename Vector<SizeType>::Type> mysort;
		typename Vector<SizeType>::Type       iperm(data_.size());
		mysort.sort(data_, iperm);
	}

	/*1. Find the largest index k such that a[k] < a[k + 1].
	 * If no such index exists, the permutation is the last permutation.
	 *
	   2. Find the largest index l such that a[k] < a[l]. Since k + 1
	   is such an index, l is well defined and satisfies k < l.

	   3. Swap a[k] with a[l].

	   4. Reverse the sequence from a[k + 1] up to and including the
	   final element a[n-1].
	   */
	bool increase()
	{
		if (data_.size() == 0)
			return false;
		int k = largestk();
		if (k < 0)
			return false;

		int l = largestl(k);
		std::swap(data_[k], data_[l]);
		SizeType c = data_.size() - 1;
		if (SizeType(k) + 1 >= data_.size() - 1)
			return true;
		typename Vector<SizeType>::Type tmp = data_;
		for (SizeType i = k + 1; i < data_.size(); i++)
			data_[c--] = tmp[i];
		return true;
	}

	SizeType operator[](SizeType i) const
	{
		// if (i>=data_.size()) throw RuntimeError("Permutations
		// error\n");
		return data_[i];
	}

	SizeType size() const { return data_.size(); }

private:

	int largestk() const
	{
		int      saved = -1;
		SizeType tot   = data_.size() - 1;
		for (SizeType i = 0; i < tot; i++) {
			if (data_[i] < data_[i + 1])
				saved = i;
		}
		return saved;
	}

	SizeType largestl(SizeType k) const
	{
		SizeType saved = 0;
		for (SizeType i = 0; i < data_.size(); i++) {
			if (data_[k] < data_[i])
				saved = i;
		}
		return saved;
	}
	typename Vector<SizeType>::Type data_;

}; // Permutations

template <typename T> std::ostream& operator<<(std::ostream& os, const Permutations<T>& ig)
{
	for (SizeType i = 0; i < ig.size(); i++)
		os << ig[i] << " ";
	return os;
}
} // namespace PsimagLite

/*@}*/
#endif // PERMUTATIONS_H_H
