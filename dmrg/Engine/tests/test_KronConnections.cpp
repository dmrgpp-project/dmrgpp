#include "KronConnections.h"
#include "KronLogger.hh"
#include <PsimagLite/CrsMatrix.h>
#include <PsimagLite/Matrix.h>
#include <catch2/catch_test_macros.hpp>
#include <complex>
#include <optional>
#include <string>
#include <vector>

namespace {

using ComplexType      = std::complex<double>;
using SparseMatrixType = PsimagLite::CrsMatrix<ComplexType>;
using VectorType       = std::vector<ComplexType>;
using BlockType        = Dmrg::MatrixDenseOrSparse<SparseMatrixType>;

struct FakeOptions {
	bool isSet(const std::string&) const { return false; }
};

struct FakeParameters {
	FakeOptions options;
	SizeType    dumperBegin = 0;
	SizeType    dumperEnd   = 0;
	std::string filename;
};

struct FakeArrayOfMatStruct {
	using MatrixDenseOrSparseType = BlockType;

	const BlockType* operator()(SizeType outPatch, SizeType inPatch) const
	{
		const auto& block = blocks[outPatch][inPatch];
		return block ? &(*block) : nullptr;
	}

	std::optional<BlockType> blocks[2][2];
};

struct FakeInitKron {
	enum WhatBasisEnum
	{
		OLD,
		NEW
	};

	using SparseMatrixType     = ::SparseMatrixType;
	using ArrayOfMatStructType = FakeArrayOfMatStruct;
	using RealType             = double;
	using ComplexOrRealType    = ComplexType;
	struct GenIjPatchType {
		using BasisType = int;
	};

	VectorType&       xout() { return xout_; }
	const VectorType& yin() const { return yin_; }
	SizeType          numberOfPatches(WhatBasisEnum) const { return 2; }
	SizeType          offsetForPatches(WhatBasisEnum, SizeType patch) const { return patch; }
	SizeType          connections() const { return 1; }
	const FakeArrayOfMatStruct& xc(SizeType) const { return x_; }
	const FakeArrayOfMatStruct& yc(SizeType) const { return y_; }
	bool                        useLowerPart() const { return false; }
	void     checks(const BlockType&, const BlockType&, SizeType, SizeType) const { }
	double   denseFlopDiscount() const { return 1.0; }
	SizeType gemmRnb() const { return 1; }
	SizeType nthreads2() const { return 1; }
	const FakeParameters& params() const { return params_; }

	FakeParameters       params_;
	FakeArrayOfMatStruct x_;
	FakeArrayOfMatStruct y_;
	VectorType           xout_;
	VectorType           yin_;
};

BlockType scalarBlock(const ComplexType& value)
{
	PsimagLite::Matrix<ComplexType> dense(1, 1);
	dense(0, 0) = value;
	return BlockType(SparseMatrixType(dense), 0.0);
}

void setBlock(FakeInitKron& initKron,
              SizeType      outPatch,
              SizeType      inPatch,
              ComplexType   left,
              ComplexType   right)
{
	initKron.x_.blocks[outPatch][inPatch].emplace(scalarBlock(left));
	initKron.y_.blocks[outPatch][inPatch].emplace(scalarBlock(right));
}

void checkVector(const VectorType& actual, const VectorType& expected)
{
	REQUIRE(actual.size() == expected.size());
	for (SizeType i = 0; i < actual.size(); ++i) {
		INFO("vector index " << i);
		CHECK(std::abs(actual[i] - expected[i]) < 1e-12);
	}
}

VectorType multiply(const ComplexType (&matrix)[2][2], const VectorType& input)
{
	VectorType output(2, ComplexType(0.0, 0.0));
	for (SizeType row = 0; row < 2; ++row)
		for (SizeType col = 0; col < 2; ++col)
			output[row] += matrix[row][col] * input[col];
	return output;
}

VectorType multiplyAdjoint(const ComplexType (&matrix)[2][2], const VectorType& input)
{
	VectorType output(2, ComplexType(0.0, 0.0));
	for (SizeType row = 0; row < 2; ++row)
		for (SizeType col = 0; col < 2; ++col)
			output[row] += std::conj(matrix[col][row]) * input[col];
	return output;
}

VectorType multiplyTranspose(const ComplexType (&matrix)[2][2], const VectorType& input)
{
	VectorType output(2, ComplexType(0.0, 0.0));
	for (SizeType row = 0; row < 2; ++row)
		for (SizeType col = 0; col < 2; ++col)
			output[row] += matrix[col][row] * input[col];
	return output;
}

} // namespace

TEST_CASE("KronConnections applies a complex non-Hermitian factor and its adjoint",
          "[KronConnections]")
{
	FakeInitKron initKron;
	setBlock(initKron, 0, 0, ComplexType(1.0, 1.0), ComplexType(2.0, -1.0));
	setBlock(initKron, 0, 1, ComplexType(2.0, 1.0), ComplexType(1.0, 2.0));
	setBlock(initKron, 1, 0, ComplexType(1.0, -2.0), ComplexType(3.0, 1.0));
	setBlock(initKron, 1, 1, ComplexType(-1.0, 1.0), ComplexType(2.0, 2.0));

	const ComplexType matrix[2][2] = { { ComplexType(3.0, 1.0), ComplexType(0.0, 5.0) },
		                           { ComplexType(5.0, -5.0), ComplexType(-4.0, 0.0) } };
	const VectorType  input { ComplexType(1.0, 1.0), ComplexType(2.0, -1.0) };

	Dmrg::KronLogger<FakeInitKron>      logger(initKron);
	Dmrg::KronConnections<FakeInitKron> forward(initKron, logger, false);
	initKron.yin_  = input;
	initKron.xout_ = { ComplexType(3.0, 1.0), ComplexType(-1.0, 2.0) };
	forward.doTask(0, 0);
	forward.doTask(1, 0);
	checkVector(initKron.xout_, { ComplexType(10.0, 15.0), ComplexType(1.0, 6.0) });

	Dmrg::KronConnections<FakeInitKron> adjoint(initKron, logger, true);
	initKron.xout_ = { ComplexType(3.0, 1.0), ComplexType(-1.0, 2.0) };
	adjoint.doTask(0, 0);
	adjoint.doTask(1, 0);
	const VectorType expectedAdjoint = multiplyAdjoint(matrix, input);
	checkVector(initKron.xout_,
	            { ComplexType(3.0, 1.0) + expectedAdjoint[0],
	              ComplexType(-1.0, 2.0) + expectedAdjoint[1] });
	const VectorType transpose = multiplyTranspose(matrix, input);
	CHECK(std::abs(expectedAdjoint[0] - transpose[0]) > 1e-12);
	CHECK(std::abs(expectedAdjoint[1] - transpose[1]) > 1e-12);
	checkVector(multiply(matrix, input), { ComplexType(7.0, 14.0), ComplexType(2.0, 4.0) });

	initKron.yin_  = input;
	initKron.xout_ = VectorType(2, ComplexType(0.0, 0.0));
	forward.doTask(0, 0);
	forward.doTask(1, 0);
	initKron.yin_  = initKron.xout_;
	initKron.xout_ = { ComplexType(3.0, 1.0), ComplexType(-1.0, 2.0) };
	adjoint.doTask(0, 0);
	adjoint.doTask(1, 0);
	const VectorType expectedFactorized = multiplyAdjoint(matrix, multiply(matrix, input));
	checkVector(initKron.xout_,
	            { ComplexType(3.0, 1.0) + expectedFactorized[0],
	              ComplexType(-1.0, 2.0) + expectedFactorized[1] });
}
