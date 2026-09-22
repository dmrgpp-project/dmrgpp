#include "MatrixMarket.hh"
#include <catch2/catch_test_macros.hpp>

#include <complex>
#include <sstream>
#include <string>

TEST_CASE("Matrix Market prints real and empty sparse matrices", "[MatrixMarket]")
{
	PsimagLite::CrsMatrix<double> matrix;
	matrix.resize(2, 3);
	matrix.setRow(0, 0);
	matrix.pushCol(2);
	matrix.pushValue(-1.5);
	matrix.setRow(1, 1);
	matrix.setRow(2, 1);

	std::ostringstream output;
	Dmrg::MatrixMarket<double>(matrix).print(output);
	CHECK(output.str()
	      == "%%MatrixMarket matrix coordinate real general\n"
	         "2 3 1\n"
	         "1 3 -1.5\n");

	PsimagLite::CrsMatrix<double> empty;
	empty.resize(2, 2);
	empty.setRow(0, 0);
	empty.setRow(1, 0);
	empty.setRow(2, 0);

	std::ostringstream emptyOutput;
	Dmrg::MatrixMarket<double>(empty).print(emptyOutput);
	CHECK(emptyOutput.str()
	      == "%%MatrixMarket matrix coordinate real general\n"
	         "2 2 0\n");
}

TEST_CASE("Matrix Market prints complex values as real-imaginary pairs", "[MatrixMarket]")
{
	using ValueType = std::complex<double>;
	PsimagLite::CrsMatrix<ValueType> matrix;
	matrix.resize(1, 1);
	matrix.setRow(0, 0);
	matrix.pushCol(0);
	matrix.pushValue(ValueType(2.0, -3.0));
	matrix.setRow(1, 1);

	std::ostringstream output;
	Dmrg::MatrixMarket<ValueType>(matrix).print(output);
	CHECK(output.str()
	      == "%%MatrixMarket matrix coordinate complex general\n"
	         "1 1 1\n"
	         "1 1 2 -3\n");
}
