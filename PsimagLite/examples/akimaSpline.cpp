// SPDX-FileCopyrightText: Copyright (c) 2009 , UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]


#include <PsimagLite/AkimaSpline.h>
#include <PsimagLite/Concurrency.h>
#include <PsimagLite/Vector.h>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

using ComplexOrRealType = double;
using VectorType        = std::vector<ComplexOrRealType>;
using AkimaSplineType   = PsimagLite::AkimaSpline<VectorType>;

void readTwoColumnData(const std::string& file, VectorType& v0, VectorType& v1)
{
	std::ifstream fin(file.c_str());
	if (!fin || !fin.good() || fin.bad())
		throw std::runtime_error("Cannot open file\n");
	while (!fin.eof()) {
		std::string s;
		fin >> s;
		if (s[0] == '#')
			continue;
		ComplexOrRealType x    = std::atof(s.c_str());
		SizeType          size = v0.size();
		if (size > 1 && x < v0[size - 1])
			break;
		v0.push_back(x);
		fin >> s;
		if (s[0] == '#')
			continue;
		v1.push_back(atof(s.c_str()));
	}
	fin.close();
}

int main(int argc, char* argv[])
{
	constexpr unsigned int nthreads = 1;
	PsimagLite::Concurrency(&argc, &argv, nthreads);

	if (argc != 5) {
		std::cerr << "USAGE: " << argv[0] << " filename start end total\n";
		return 1;
	}

	VectorType x, s;
	readTwoColumnData(argv[1], x, s);

	AkimaSplineType   akimaSpline(x, s);
	ComplexOrRealType xstart = std::atof(argv[2]);
	ComplexOrRealType xend   = std::atof(argv[3]);
	SizeType          total  = std::atoi(argv[4]);
	ComplexOrRealType xstep  = (xend - xstart) / total;

	for (ComplexOrRealType x = xstart; x < xend; x += xstep) {
		std::cout << x << " " << akimaSpline(x) << "\n";
	}
}
