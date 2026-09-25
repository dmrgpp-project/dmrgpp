// SPDX-FileCopyrightText: Copyright (c) 2009 , UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup PsimagLite */
/*@{*/

/*! \file Stack.h
 *
 *  std::stack companions
 */

#ifndef PSIMAGLITE_STACK_H_
#define PSIMAGLITE_STACK_H_
#include "AllocatorCpu.h"
#include "Vector.h"
#include <stack>

namespace PsimagLite {

template <typename T> class Stack {

	using DequeType_ = std::deque<T, typename Allocator<T>::Type>;

public:

	using Type = std::stack<T, DequeType_>;
}; // class Stack

template <typename T> class IsStackLike {
public:

	enum
	{
		True = false
	};
};

template <typename T> class IsStackLike<std::stack<T, std::deque<T, typename Allocator<T>::Type>>> {
public:

	enum
	{
		True = true
	};
};

template <typename StackType>
typename EnableIf<IsStackLike<StackType>::True, std::ostream>::Type& operator<<(std::ostream&    os,
                                                                                const StackType& st)
{
	StackType st2 = st;
	os << st2.size() << "\n";
	while (!st2.empty()) {
		typename StackType::value_type x = st2.top();
		os << x << "\n";
		st2.pop();
	}
	return os;
}

template <typename StackType>
typename EnableIf<IsStackLike<StackType>::True, std::istream>::Type& operator>>(std::istream& is,
                                                                                StackType&    x)
{
	using ValueType = typename StackType::value_type;
	typename Vector<ValueType>::Type tmpVec;
	is >> tmpVec;
	for (int i = tmpVec.size() - 1; i >= 0; i--) {
		x.push(tmpVec[i]);
	}
	return is;
}

} // namespace PsimagLite

/*@}*/
#endif // PSIMAGLITE_STACK_H_
