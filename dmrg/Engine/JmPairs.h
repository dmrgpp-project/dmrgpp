// SPDX-FileCopyrightText: Copyright (c) 2009-2017, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 4.]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file JmPairs.h
 *
 *  This is a "vector" of (2j,m+j) pairs.
 *  Repeated entries are stored only once
 *  Provides a transparent access as if it were a normal vector
 *
 */
#ifndef JMPAIRS_HEADER_H
#define JMPAIRS_HEADER_H

#include "Utils.h"
#include <PsimagLite/PsimagLite.h>
#include <algorithm>

namespace Dmrg {

template <typename PairType_> class JmPairs {

	using VectorSizeType = typename PsimagLite::Vector<SizeType>::Type;

public:

	using PairType       = PairType_;
	using value_type     = PairType;
	using VectorPairType = typename PsimagLite::Vector<PairType>::Type;

	//! indices_[alpha] = jm
	void push(const PairType& jm, SizeType)
	{
		int x = PsimagLite::indexOrMinusOne(jmPairs_, jm);

		if (x < 0) {
			jmPairs_.push_back(jm);
			x = jmPairs_.size() - 1;
		}

		indices_.push_back(x);
	}

	const PairType& operator[](SizeType alpha) const
	{
		assert(alpha < indices_.size());
		assert(indices_[alpha] < jmPairs_.size());
		return jmPairs_[indices_[alpha]];
	}

	JmPairs<PairType>& operator=(const VectorPairType& jmpairs)
	{
		SizeType n = jmpairs.size();
		jmPairs_.clear();
		indices_.resize(n);
		for (SizeType i = 0; i < n; ++i) {
			int x = PsimagLite::indexOrMinusOne(jmPairs_, jmpairs[i]);
			if (x < 0) {
				jmPairs_.push_back(jmpairs[i]);
				x = jmPairs_.size() - 1;
			}

			indices_[i] = x;
		}

		return *this;
	}

	void clear()
	{
		jmPairs_.clear();
		indices_.clear();
	}

	void reorder(const VectorSizeType& permutation) { utils::reorder(indices_, permutation); }

	void truncate(const VectorSizeType& removedIndices)
	{
		utils::truncateVector(indices_, removedIndices);
		VectorSizeType unusedPairs;
		findUnusedJmPairs(unusedPairs);
		removeUnusedPairs(unusedPairs);
	}

	template <typename Op> void maxFirst(SizeType& maxvalue)
	{
		Op f;
		for (SizeType i = 0; i < jmPairs_.size(); i++) {
			if (f(jmPairs_[i].first, maxvalue)) {
				maxvalue = jmPairs_[i].first;
			}
		}
	}

	SizeType size() const { return indices_.size(); }

	template <typename IoOutputter>
	void
	write(IoOutputter&       io,
	      PsimagLite::String prefix,
	      typename PsimagLite::EnableIf<PsimagLite::IsOutputLike<IoOutputter>::True, int>::Type
	      = 0) const
	{
		io.write(jmPairs_, prefix + "su2JmPairs");
		io.write(indices_, prefix + "su2JmIndices");
	}

	template <typename IoInputter>
	void
	read(IoInputter&        io,
	     PsimagLite::String prefix,
	     typename PsimagLite::EnableIf<PsimagLite::IsInputLike<IoInputter>::True, int>::Type
	     = 0)
	{
		io.read(jmPairs_, prefix + "su2JmPairs");
		io.read(indices_, prefix + "su2JmIndices");
	}

	friend std::ostream& operator<<(std::ostream& os, JmPairs<PairType> jmPairs)
	{
		for (SizeType i = 0; i < jmPairs.size(); i++)
			os << "jmPair[" << i << "]=" << jmPairs[i] << "\n";
		return os;
	}

private:

	void findUnusedJmPairs(VectorSizeType& unusedPairs)
	{
		for (SizeType i = 0; i < jmPairs_.size(); i++)
			if (isUnusedPair(i))
				unusedPairs.push_back(i);
	}

	void removeUnusedPairs(const VectorSizeType& unusedPairs)
	{
		SizeType       counter = 0;
		VectorSizeType neworder(jmPairs_.size());
		VectorPairType tmpVector(jmPairs_.size() - unusedPairs.size());

		for (SizeType i = 0; i < jmPairs_.size(); i++) {
			if (PsimagLite::indexOrMinusOne(unusedPairs, i) >= 0)
				continue;
			tmpVector[counter] = jmPairs_[i];
			neworder[i]        = counter;
			counter++;
		}

		jmPairs_ = tmpVector;
		VectorSizeType tmpVector2(indices_.size());
		for (SizeType i = 0; i < indices_.size(); i++)
			tmpVector2[i] = neworder[indices_[i]];
		indices_ = tmpVector2;
	}

	bool isUnusedPair(SizeType ind)
	{
		for (SizeType i = 0; i < indices_.size(); i++)
			if (indices_[i] == ind)
				return false;
		return true;
	}

	VectorPairType jmPairs_;
	VectorSizeType indices_;
}; // JmPairs
} // namespace Dmrg
/*@}*/
#endif
