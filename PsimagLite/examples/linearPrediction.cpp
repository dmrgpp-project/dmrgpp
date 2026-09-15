// SPDX-FileCopyrightText: Copyright (c) 2009 , UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]


#include <PsimagLite/Io/IoSimple.h>
#include <PsimagLite/LinearPrediction.h>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

typedef double                                  FieldType;
typedef PsimagLite::LinearPrediction<FieldType> LinearPredictionType;

void usage(const char* progName)
{
	std::cerr << "Usage: " << progName << " -f file -l label -p n -q q[>1]\n";
}

int main(int argc, char* argv[])
{
	int                opt   = 0;
	PsimagLite::String file  = "";
	PsimagLite::String label = "";
	SizeType           p     = 0;
	SizeType           q     = 2;
	while ((opt = getopt(argc, argv, "f:p:l:q:")) != -1) {
		switch (opt) {
		case 'f':
			file = optarg;
			break;
		case 'l':
			label = optarg;
			break;
		case 'p':
			p = atoi(optarg);
			break;
		case 'q':
			q = atoi(optarg);
			break;
		default:
			usage(argv[0]);
			return 1;
		}
	}
	// sanity checks:
	if (file == "" || label == "" || p == 0 || q < 2) {
		usage(argv[0]);
		return 1;
	}

	PsimagLite::IoSimple::In            io(file);
	PsimagLite::Vector<FieldType>::Type y;
	io.read(y, label);
	SizeType n = y.size();
	std::cout << "#Found " << n << " points in file " << file << "\n";
	LinearPredictionType linearPrediction(y, q);
	linearPrediction.predict(q);

	for (SizeType i = 0; i < p; i++) {
		linearPrediction.linearPredictionfunction(y, q);
		linearPrediction.predict(q);
	}

	for (SizeType i = 0; i < p + n; i++) {
		std::cout << linearPrediction(i) << "\n";
	}
}

// Test Function in series.txt is: f[i] =
// 0.5*(cos((i-1)*pi/100)+cos((i-1)*pi/20))*exp(-(i-1)/100)
