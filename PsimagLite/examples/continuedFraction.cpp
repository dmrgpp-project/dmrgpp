// SPDX-FileCopyrightText: Copyright (c) 2009 , UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]


#include <PsimagLite/ContinuedFraction.h>
#include <PsimagLite/Io/IoSimple.h>
#include <PsimagLite/TridiagonalMatrix.h>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

using namespace PsimagLite;
typedef double                      RealType;
typedef ContinuedFraction<RealType> ContinuedFractionType;

void usage(const char* progName)
{
	std::cerr << "Usage: " << progName << " -f file  -b omega1";
	std::cerr << " -t total -s omegaStep -d delta\n";
	std::cerr << "Conditions: total > 0 omegaStep>0 delta>0\n";
}

int main(int argc, char* argv[])
{
	int      opt    = 0;
	String   file   = "";
	RealType wbegin = 0;
	RealType wstep  = 0;
	RealType delta  = 0;
	SizeType total  = 0;

	while ((opt = getopt(argc, argv, "f:b:t:s:d:")) != -1) {
		switch (opt) {
		case 'f':
			file = optarg;
			break;
		case 'b':
			wbegin = atof(optarg);
			break;
		case 't':
			total = atoi(optarg);
			break;
		case 's':
			wstep = atof(optarg);
			break;
		case 'd':
			delta = atof(optarg);
			break;
		default:
			usage(argv[0]);
			return 1;
		}
	}
	// sanity checks:
	if (file == "" || total <= 0 || wstep <= 0 || delta <= 0) {
		usage(argv[0]);
		return 1;
	}

	IoSimple::In                                 io(file);
	ContinuedFractionType                        cf(io);
	typename ContinuedFractionType::PlotDataType v;
	RealFrequencyRange<double>                   plotParams(wbegin, wstep, total, delta);
	cf.plot(v, plotParams);
	for (SizeType x = 0; x < v.size(); x++) {
		std::cout << v[x].first << " " << PsimagLite::real(v[x].second);
		std::cout << " " << PsimagLite::imag(v[x].second) << "\n";
	}
}
