// SPDX-FileCopyrightText: Copyright (c) 2026, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later

#include <PsimagLite/CrsMatrix.h>
#include <PsimagLite/Matrix.h>
#include <catch2/catch_test_macros.hpp>

#include <complex>
#include <random>
#include <vector>

namespace {

template <typename T>
PsimagLite::CrsMatrix<T> makeCrs(SizeType rows, SizeType cols, const std::vector<T>& values)
{
	REQUIRE(values.size() == rows * cols);
	PsimagLite::Matrix<T> matrix(rows, cols);
	for (SizeType i = 0; i < rows; ++i)
		for (SizeType j = 0; j < cols; ++j)
			matrix(i, j) = values[i * cols + j];

	return PsimagLite::CrsMatrix<T>(matrix);
}

template <typename T>
PsimagLite::Matrix<T> denseMultiply(const PsimagLite::Matrix<T>& a, const PsimagLite::Matrix<T>& b)
{
	REQUIRE(a.cols() == b.rows());
	PsimagLite::Matrix<T> result(a.rows(), b.cols());
	for (SizeType i = 0; i < a.rows(); ++i)
		for (SizeType j = 0; j < b.cols(); ++j)
			for (SizeType k = 0; k < a.cols(); ++k)
				result(i, j) += a(i, k) * b(k, j);

	return result;
}

template <typename T>
void requireEqual(const PsimagLite::CrsMatrix<T>& actual, const PsimagLite::Matrix<T>& expected)
{
	REQUIRE(actual.rows() == expected.rows());
	REQUIRE(actual.cols() == expected.cols());
	const PsimagLite::Matrix<T> dense = actual.toDense();
	for (SizeType i = 0; i < expected.rows(); ++i)
		for (SizeType j = 0; j < expected.cols(); ++j)
			REQUIRE(dense(i, j) == expected(i, j));
}

} // namespace

TEST_CASE("CrsMatrix multiply handles rectangular matrices and sorts output columns",
          "[CrsMatrix][multiply]")
{
	const auto a = makeCrs<double>(2, 3, { 1, 0, 2, 0, 3, 4 });
	const auto b = makeCrs<double>(3, 4, { 0, 5, 0, 1, 7, 0, 8, 0, 0, 9, 10, 0 });

	PsimagLite::CrsMatrix<double> result;
	PsimagLite::multiply(result, a, b);

	requireEqual(result, makeCrs<double>(2, 4, { 0, 23, 20, 1, 21, 36, 64, 0 }).toDense());
	for (SizeType i = 0; i < result.rows(); ++i)
		for (int k = result.getRowPtr(i) + 1; k < result.getRowPtr(i + 1); ++k)
			REQUIRE(result.getCol(k - 1) < result.getCol(k));
}

TEST_CASE("CrsMatrix multiply removes entries cancelled to exact zero", "[CrsMatrix][multiply]")
{
	const auto a = makeCrs<double>(1, 2, { 1, 1 });
	const auto b = makeCrs<double>(2, 1, { 1, -1 });

	PsimagLite::CrsMatrix<double> result;
	PsimagLite::multiply(result, a, b);

	REQUIRE(result.rows() == 1);
	REQUIRE(result.cols() == 1);
	REQUIRE(result.nonZeros() == 0);
}

TEST_CASE("CrsMatrix multiply preserves empty rows and empty products", "[CrsMatrix][multiply]")
{
	const auto a = makeCrs<double>(3, 2, { 1, 0, 0, 0, 0, 2 });
	const auto b = makeCrs<double>(2, 3, { 0, 0, 0, 3, 0, 0 });

	PsimagLite::CrsMatrix<double> result;
	PsimagLite::multiply(result, a, b);

	requireEqual(result, makeCrs<double>(3, 3, { 0, 0, 0, 0, 0, 0, 6, 0, 0 }).toDense());
	REQUIRE(result.getRowPtr(0) == result.getRowPtr(1));
	REQUIRE(result.getRowPtr(1) == result.getRowPtr(2));
}

TEST_CASE("CrsMatrix multiply preserves a matrix multiplied by the identity",
          "[CrsMatrix][multiply]")
{
	const auto matrix   = makeCrs<double>(3, 3, { 0, 2, 0, -1, 0, 4, 5, 0, 0 });
	const auto identity = makeCrs<double>(3, 3, { 1, 0, 0, 0, 1, 0, 0, 0, 1 });

	PsimagLite::CrsMatrix<double> result;
	PsimagLite::multiply(result, matrix, identity);

	REQUIRE(result == matrix);
}

TEST_CASE("CrsMatrix multiply supports mixed input scalar types", "[CrsMatrix][multiply]")
{
	using Complex = std::complex<double>;
	const auto a  = makeCrs<double>(2, 2, { 1, 2, 0, -1 });
	const auto b  = makeCrs<Complex>(
            2, 2, { Complex(0, 1), Complex(3, 0), Complex(2, -1), Complex(0, 4) });

	PsimagLite::CrsMatrix<Complex> result;
	PsimagLite::multiply(result, a, b);

	requireEqual(result,
	             makeCrs<Complex>(
	                 2, 2, { Complex(4, -1), Complex(3, 8), Complex(-2, 1), Complex(0, -4) })
	                 .toDense());
}

TEST_CASE("CrsMatrix multiply rejects incompatible dimensions", "[CrsMatrix][multiply]")
{
	const auto                    a = makeCrs<double>(1, 2, { 0, 0 });
	const auto                    b = makeCrs<double>(3, 1, { 0, 0, 0 });
	PsimagLite::CrsMatrix<double> result;

	REQUIRE_THROWS(PsimagLite::multiply(result, a, b));
}

TEST_CASE("CrsMatrix multiply permits output to alias either input", "[CrsMatrix][multiply]")
{
	const auto originalA = makeCrs<double>(2, 2, { 1, 2, 0, 3 });
	const auto originalB = makeCrs<double>(2, 2, { 4, 0, 5, 6 });
	const auto expected  = denseMultiply(originalA.toDense(), originalB.toDense());

	SECTION("Output aliases A")
	{
		auto a = originalA;
		PsimagLite::multiply(a, a, originalB);
		requireEqual(a, expected);
	}

	SECTION("Output aliases B")
	{
		auto b = originalB;
		PsimagLite::multiply(b, originalA, b);
		requireEqual(b, expected);
	}
}

TEST_CASE("CrsMatrix multiply agrees with dense multiplication on random small matrices",
          "[CrsMatrix][multiply]")
{
	std::mt19937                       rng(1729);
	std::uniform_int_distribution<int> value(-2, 2);
	std::bernoulli_distribution        keep(0.4);

	for (SizeType iteration = 0; iteration < 40; ++iteration) {
		PsimagLite::Matrix<double> denseA(4, 5);
		PsimagLite::Matrix<double> denseB(5, 3);
		for (SizeType i = 0; i < denseA.rows(); ++i)
			for (SizeType j = 0; j < denseA.cols(); ++j)
				denseA(i, j) = keep(rng) ? value(rng) : 0;
		for (SizeType i = 0; i < denseB.rows(); ++i)
			for (SizeType j = 0; j < denseB.cols(); ++j)
				denseB(i, j) = keep(rng) ? value(rng) : 0;

		const PsimagLite::CrsMatrix<double> a(denseA);
		const PsimagLite::CrsMatrix<double> b(denseB);
		PsimagLite::CrsMatrix<double>       result;
		PsimagLite::multiply(result, a, b);

		requireEqual(result, denseMultiply(denseA, denseB));
	}
}
