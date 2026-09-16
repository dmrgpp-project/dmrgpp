// SPDX-FileCopyrightText: Copyright (c) 2009-2016-2018, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [by G.A., Oak Ridge National Laboratory]

#include "Utils.h"

// Utility functions that are still needed
namespace utils {

PsimagLite::String pathPrepend(PsimagLite::String pre, PsimagLite::String pathname)
{
	bool               addDotDot = false;
	PsimagLite::String path1("");
	if (pathname.length() > 2 && pathname[0] == '.' && pathname[1] == '.') {
		addDotDot = true;
		path1     = pathname.substr(2, pathname.length());
	}

	if (pathname.length() > 1 && pathname[0] == '/')
		path1 = pathname;

	if (path1 == "")
		return pre + pathname;

	size_t index = path1.find_last_of("/");

	index++;
	PsimagLite::String ret = path1.substr(0, index) + pre + path1.substr(index, path1.length());
	return (addDotDot) ? ".." + ret : ret;
}

SizeType exactDivision(SizeType a, SizeType b)
{
	SizeType c = static_cast<SizeType>(a / b);
	if (c * b != a)
		throw PsimagLite::RuntimeError("exactDivision expected\n");

	return c;
}

SizeType bitSizeOfInteger(SizeType x)
{
	SizeType counter = 0;
	while (x) {
		counter++;
		x >>= 1;
	}

	return counter;
}

SizeType powUint(SizeType x, SizeType n)
{
	if (n == 0)
		return 1;
	if (n == 1)
		return x;
	SizeType ret = x;
	for (SizeType i = 1; i < n; ++i)
		ret *= x;
	return ret;
}

} // namespace utils
