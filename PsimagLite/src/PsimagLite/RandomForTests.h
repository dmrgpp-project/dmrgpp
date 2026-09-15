// SPDX-FileCopyrightText: Copyright (c) 2009-2013, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]


/** \ingroup PsimagLite */
/*@{*/

/*! \file RandomForTests.h
 *
 *  THIS RNG IS FOR TESTS, DO NOT USE FOR PRODUCTION, IT'S NOT RANDOM ENOUGH!!
 *
 */

#ifndef RANDOM_FOR_TESTS_H
#define RANDOM_FOR_TESTS_H
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

namespace PsimagLite {
template <typename T> class RandomForTests {

public:

	using LongType = long int;
	typedef T value_type; // legacy name

	RandomForTests(int seed)
	    : next_(seed)
	{ }

	T operator()()
	{
		next_ = 16807 * (next_ % 127773) - 2836 * (next_ / 127773);
		if (next_ <= 0)
			next_ += 2147483647;
		if (next_ >= 2147483647)
			next_ = 1;
		return static_cast<T>(next_) / 2147483647.0;
	}

private:

	int next_;
}; // RandomForTests
} // namespace PsimagLite

/*@}*/
#endif // RANDOM_FOR_TESTS_H
