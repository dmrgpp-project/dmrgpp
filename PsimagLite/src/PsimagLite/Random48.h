// SPDX-FileCopyrightText: Copyright (c) 2009-2013, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]


#ifndef RANDOM48_H
#define RANDOM48_H
#include "Vector.h"
#include <cstdlib>
#include <iostream>

namespace PsimagLite {
template <typename T> class Random48 {

public:

	using LongType = long int;
	typedef T value_type; // legacy name

	Random48(LongType seed, SizeType rank = 0, SizeType nprocs = 1)
	{
		srand48(seed);
		Vector<LongType>::Type vOfSeeds(nprocs);
		for (SizeType i = 0; i < vOfSeeds.size(); i++)
			vOfSeeds[i] = static_cast<LongType>(10000.0 * drand48());
		seed_ = vOfSeeds[rank];
		srand48(seed_);
	}

	T random() const // deprecated!!! use operator() instead
	{
		return static_cast<T>(drand48());
	}

	T operator()() const { return static_cast<T>(drand48()); }

	LongType seed() const { return seed_; }

	void seed(const LongType& seed) { seed_ = seed; }

private:

	LongType seed_;
}; // Random48
} // namespace PsimagLite

#endif // RANDOM48_H
