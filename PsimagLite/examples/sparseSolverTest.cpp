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

#include <PsimagLite/CrsMatrix.h>
#include <PsimagLite/LanczosSolver.h>
#include <PsimagLite/ParametersForSolver.h>
#include <PsimagLite/PsimagLite.h>
#include <PsimagLite/Random48.h>

using namespace PsimagLite;

using RealType          = double;
using ComplexOrRealType = double;

using SparseMatrixType  = CrsMatrix<ComplexOrRealType>;
using MatrixSolverType  = MatrixSolverBase<SparseMatrixType>;
using LanczosSolverType = LanczosSolver<SparseMatrixType>;

void usage(const char* progName)
{
	std::cerr << "Usage: " << progName
	          << " -n rank [-x] [-d ] [-c max_columns] [-m max_value] [-r seed]\n";
	exit(1);
}

int main(int argc, char* argv[])
{
	constexpr unsigned int  nthreads = 1;
	PsimagLite::Concurrency concurrency(&argc, &argv, nthreads);

	int      opt         = 0;
	SizeType n           = 0;
	RealType max_value   = 0;
	SizeType max_col     = 0;
	SizeType seed        = 0;
	bool     lota_memory = true;

	while ((opt = getopt(argc, argv, "n:c:m:r:x")) != -1) {
		switch (opt) {
		case 'n':
			n = atoi(optarg);
			break;
		case 'c':
			max_col = atoi(optarg);
			break;
		case 'm':
			max_value = atof(optarg);
			break;
		case 'r':
			seed = atoi(optarg);
			break;
		case 'x':
			lota_memory = false;
			break;
		default:
			usage(argv[0]);
			return 1;
		}
	}

	// sanity checks
	if (n == 0)
		usage(argv[0]);
	if (max_col == 0)
		max_col = 1 + SizeType(0.1 * n);
	if (PsimagLite::norm(max_value) < 1e-6)
		max_value = 1.0;
	if (seed == 0)
		seed = 3443331;

	// create a random matrix:
	Random48<RealType> random(seed);
	SparseMatrixType   sparse(n, n);
	Vector<bool>::Type seen_this_column(n);
	SizeType           counter = 0;
	for (SizeType i = 0; i < n; i++) {
		sparse.setRow(i, counter);
		// random vector:
		SizeType x = 1 + SizeType(random() * max_col);
		for (SizeType j = 0; j < seen_this_column.size(); j++)
			seen_this_column[j] = false;
		for (SizeType j = 0; j < x; j++) {
			SizeType col = SizeType(random() * n);
			if (seen_this_column[col])
				continue;
			seen_this_column[col] = true;
			ComplexOrRealType val = random() * max_value;

			sparse.pushValue(val);
			sparse.pushCol(col);
			counter++;
		}
	}
	sparse.setRow(n, counter);
	sparse.checkValidity();
	// symmetrize:
	SparseMatrixType sparse2;
	transposeConjugate(sparse2, sparse);
	sparse += sparse2;
	sparse.checkValidity();
	assert(isHermitian(sparse));

	// sparse solver setup
	PsimagLite::ParametersForSolver<RealType> params;
	params.lotaMemory = lota_memory;
	LanczosSolverType lanczos_solver(sparse, params);

	// diagonalize matrix
	RealType                       gs_energy = 0;
	std::vector<ComplexOrRealType> gs_vector(n);

	std::vector<ComplexOrRealType> initial(n);
	PsimagLite::fillRandom(initial);
	lanczos_solver.computeOneState(gs_energy, gsVector, initial, 0);

	std::cout << "Energy=" << gs_energy << "\n";
}
