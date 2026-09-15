// SPDX-FileCopyrightText: Copyright (c) 2009 , UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]


#include <PsimagLite/ContinuedFraction.h>
#include <PsimagLite/ContinuedFractionCollection.h>
#include <PsimagLite/Io/IoSimple.h>
#include <PsimagLite/Matsubaras.h>
#include <PsimagLite/RealFrequencyRange.hh>
#include <PsimagLite/TridiagonalMatrix.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

using namespace PsimagLite;
typedef double                                RealType;
typedef ContinuedFraction<RealType>           ContinuedFractionType;
typedef ContinuedFractionCollection<RealType> ContinuedFractionCollectionType;

void usage(const char* progName)
{
	std::cerr << "Usage: " << progName << " -f file  -b omega1";
	std::cerr << " -t total -s omegaStep -d delta -B beta -m matsubaras\n";
}

template <typename PlotParamsType>
void plotAll(const ContinuedFractionCollectionType& cfCollection, const PlotParamsType& params)
{
	ContinuedFractionCollectionType::PlotDataType v;
	cfCollection.plot(v, params);
	std::cout.precision(12);
	for (SizeType x = 0; x < v.size(); x++) {
		std::cout << v[x].first << " " << PsimagLite::imag(v[x].second);
		std::cout << " " << PsimagLite::real(v[x].second) << "\n";
	}
}

template <typename PlotParamsType>
void plotOneByOne(const ContinuedFractionCollectionType& cfCollection, const PlotParamsType& params)
{

	std::cout.precision(12);
	for (SizeType i = 0; i < cfCollection.size(); i++) {
		ContinuedFractionCollectionType::PlotDataType v;
		cfCollection.plotOne(i, v, params);
		for (SizeType x = 0; x < v.size(); x++) {
			std::cout << v[x].first << " " << PsimagLite::imag(v[x].second);
			std::cout << " " << PsimagLite::real(v[x].second) << "\n";
		}
	}
}

template <typename SomeParamsType>
void plot(const SomeParamsType&                  params,
          const ContinuedFractionCollectionType& cfCollection,
          bool                                   oneByOne)
{
	if (!oneByOne)
		plotAll(cfCollection, params);
	else
		plotOneByOne(cfCollection, params);
}

int main(int argc, char* argv[])
{
	int      opt      = 0;
	String   file     = "";
	RealType wbegin   = 0;
	SizeType total    = 0;
	RealType wstep    = 0;
	RealType delta    = 0;
	RealType beta     = 0.0;
	bool     oneByOne = false;
	while ((opt = getopt(argc, argv, "f:b:t:s:d:B:1")) != -1) {
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
		case '1':
			oneByOne = true;
			break;
		case 'B':
			beta = atof(optarg);
			break;
		default:
			usage(argv[0]);
			return 1;
		}
	}
	// sanity checks:
	bool not_valid_real = (beta > 0 && (wstep <= 0 || delta <= 0));

	if (file == "" || not_valid_real || beta < 0 || total == 0) {
		usage(argv[0]);
		return 1;
	}

	IoSimple::In                    io(file);
	ContinuedFractionCollectionType cfCollection(io);

	bool is_matsubaras = (beta > 0);
	if (is_matsubaras) {
		PsimagLite::Matsubaras<RealType> matsubaras(beta, total, delta);
		plot(matsubaras, cfCollection, oneByOne);
	} else {
		PsimagLite::RealFrequencyRange<RealType> real_freq(wbegin, wstep, total, delta);
		plot<>(real_freq, cfCollection, oneByOne);
	}
}
