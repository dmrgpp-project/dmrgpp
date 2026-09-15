// SPDX-FileCopyrightText: Copyright (c) 2009-2016-2018, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [by G.A., Oak Ridge National Laboratory]

#ifndef UTILS_HEADER_H
#define UTILS_HEADER_H

#include <PsimagLite/CrsMatrix.h>
#include <PsimagLite/Vector.h>

namespace std {

template <class T1, class T2> ostream& operator<<(std::ostream& os, const pair<T1, T2>& p)
{
	os << p.first << " " << p.second;
	return os;
}
} // namespace std

// Utility functions that are still needed
namespace utils {

PsimagLite::String pathPrepend(PsimagLite::String, PsimagLite::String);

template <template <typename, typename> class SomeVectorTemplate,
          typename SomeAllocator1Type,
          typename SomeAllocator2Type,
          typename T>
typename PsimagLite::
    EnableIf<PsimagLite::IsVectorLike<SomeVectorTemplate<T, SomeAllocator1Type>>::True, void>::Type
    reorder(SomeVectorTemplate<T, SomeAllocator1Type>&              v,
            const SomeVectorTemplate<SizeType, SomeAllocator2Type>& permutation)
{
	SomeVectorTemplate<T, SomeAllocator1Type> tmpVector(v.size());
	for (SizeType i = 0; i < v.size(); i++)
		tmpVector[i] = v[permutation[i]];
	v = tmpVector;
}

template <typename SomeType>
void reorder(PsimagLite::Matrix<SomeType>& v, const PsimagLite::Vector<SizeType>::Type& permutation)
{
	PsimagLite::Matrix<SomeType> tmpVector(v.n_row(), v.n_col());
	for (SizeType i = 0; i < v.n_row(); i++)
		for (SizeType j = 0; j < v.n_col(); j++)
			tmpVector(i, j) = v(permutation[i], permutation[j]);
	v = tmpVector;
}

//! A = B union C
template <typename Block> void blockUnion(Block& A, Block const& B, Block const& C)
{
	A = B;
	for (SizeType i = 0; i < C.size(); i++)
		A.push_back(C[i]);
}

template <typename SomeVectorType>
typename PsimagLite::EnableIf<PsimagLite::IsVectorLike<SomeVectorType>::True, void>::Type
truncateVector(SomeVectorType& v, const PsimagLite::Vector<SizeType>::Type& removedIndices)
{
	SomeVectorType tmpVector;
	for (SizeType i = 0; i < v.size(); i++) {
		if (PsimagLite::indexOrMinusOne(removedIndices, i) >= 0)
			continue;
		tmpVector.push_back(v[i]);
	}

	v = tmpVector;
}

template <typename SomeVectorType>
static typename PsimagLite::EnableIf<PsimagLite::IsVectorLike<SomeVectorType>::True, void>::Type
fillFermionicSigns(SomeVectorType&                       fermionicSigns,
                   const PsimagLite::Vector<bool>::Type& signs,
                   int                                   f)
{
	using ValueType = typename SomeVectorType::value_type;
	fermionicSigns.resize(signs.size());
	for (SizeType i = 0; i < fermionicSigns.size(); i++)
		fermionicSigns[i] = (signs[i]) ? static_cast<ValueType>(f) : 1;
}

SizeType exactDivision(SizeType, SizeType);

SizeType bitSizeOfInteger(SizeType);

SizeType powUint(SizeType, SizeType);

} // namespace utils
/*@}*/
#endif
