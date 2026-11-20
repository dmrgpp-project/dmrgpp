// SPDX-FileCopyrightText: Copyright (c) 2009 , UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]

#include <PsimagLite/ContinuedFraction.h>
#include <PsimagLite/ContinuedFractionCollection.h>
#include <PsimagLite/Io/IoSimple.h>
#include <PsimagLite/TridiagonalMatrix.h>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

using namespace PsimagLite;
typedef double                                RealType;
typedef ContinuedFraction<RealType>           ContinuedFractionType;
typedef ContinuedFractionCollection<RealType> ContinuedFractionCollectionType;

void usage(const char* progName) { std::cerr << "Usage: " << progName << " file1 file2\n"; }

int main(int argc, char* argv[])
{
	if (argc < 2) {
		usage(argv[0]);
		return 1;
	}

	ContinuedFractionCollectionType cfCollection(PsimagLite::FreqEnum::REAL);

	String s = "#Avector";
	for (int x = 1; x < argc; x++) {
		IoSimple::In io(argv[x]);
		io.advance(s, IoSimple::In::LAST_INSTANCE);
		ContinuedFractionType cf(io);
		cfCollection.push(cf);
	}

	IoSimple::Out ioOut(std::cout);
	ioOut.setPrecision(12);
	cfCollection.write(ioOut);
}
