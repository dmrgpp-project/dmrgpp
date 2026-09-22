// BEGIN LICENSE BLOCK
/*
Copyright (c) 2009 , UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 1.0.0]

*********************************************************
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED.

Please see full open source license included in file LICENSE.
*********************************************************

*/
// END LICENSE BLOCK

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

	ContinuedFractionCollectionType cf_collection(PsimagLite::FreqEnum::REAL);

	String s = "#Avector";
	for (int x = 1; x < argc; x++) {
		IoSimple::In io(argv[x]);
		io.advance(s, IoSimple::In::LAST_INSTANCE);
		ContinuedFractionType cf(io);
		cf_collection.push(cf);
	}

	IoSimple::Out io_out(std::cout);
	io_out.setPrecision(12);
	cf_collection.write(io_out);
}
