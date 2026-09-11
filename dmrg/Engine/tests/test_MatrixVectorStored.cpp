#include "MatrixVectorStored.h"
#include <PsimagLite/CrsMatrix.h>
#include <PsimagLite/Matrix.h>
#include <catch2/catch_test_macros.hpp>
#include <complex>
#include <string>
#include <vector>

namespace {

using ComplexType      = std::complex<double>;
using SparseMatrixType = PsimagLite::CrsMatrix<ComplexType>;

struct FakeOptions {
	using Type = FakeOptions;

	bool isSet(const std::string& option) const { return option == "LdaggerL" && ldaggerL; }

	bool ldaggerL = false;
};

struct FakeParameters {
	using OptionsType = FakeOptions;

	FakeOptions options;
	int         maxMatrixRankStored = 10;
};

struct FakeAux {
	SizeType m() const { return 0; }
};

struct FakeSuper {
	SizeType partition(SizeType index) const { return 2 * index; }
};

struct FakeLeftRightSuper {
	const FakeSuper& super() const { return super_; }

	FakeSuper super_;
};

struct FakeModelHelper {
	using SparseMatrixType = ::SparseMatrixType;
	using RealType         = double;
	using Aux              = FakeAux;

	const FakeLeftRightSuper& leftRightSuper() const { return leftRightSuper_; }

	FakeLeftRightSuper leftRightSuper_;
};

struct FakeHamiltonianConnection {
	FakeHamiltonianConnection(const SparseMatrixType& matrix)
	    : matrix_(matrix)
	{ }

	void fullHamiltonian(SparseMatrixType& matrix, const FakeAux&, bool) const
	{
		matrix = matrix_;
	}

	const FakeModelHelper& modelHelper() const { return modelHelper_; }

	SparseMatrixType matrix_;
	FakeModelHelper  modelHelper_;
};

struct FakeModel {
	using HamiltonianConnectionType = FakeHamiltonianConnection;
	using ParametersType            = FakeParameters;
	using ModelHelperType           = FakeModelHelper;

	const FakeParameters& params() const { return params_; }
	bool                  isHermitian() const { return false; }

	FakeParameters params_;
};

struct FakeTypes {
	using ModelType        = FakeModel;
	using ModelHelperType  = FakeModelHelper;
	using RealType         = double;
	using SparseMatrixType = ::SparseMatrixType;
};

SparseMatrixType makeMatrix()
{
	PsimagLite::Matrix<ComplexType> matrix(2, 2);
	matrix(0, 0) = ComplexType(1.0, 1.0);
	matrix(0, 1) = ComplexType(2.0, 0.0);
	matrix(1, 0) = ComplexType(0.0, 3.0);
	matrix(1, 1) = ComplexType(-1.0, 0.0);
	return SparseMatrixType(matrix);
}

void checkVector(const std::vector<ComplexType>& actual, const std::vector<ComplexType>& expected)
{
	REQUIRE(actual.size() == expected.size());
	for (SizeType i = 0; i < actual.size(); ++i) {
		INFO("vector index " << i);
		CHECK(std::abs(actual[i] - expected[i]) < 1e-12);
	}
}

} // namespace

TEST_CASE("MatrixVectorStored accepts a test type bundle", "[MatrixVectorStored]")
{
	FakeModel                 model;
	FakeHamiltonianConnection connection(makeMatrix());
	FakeAux                   aux;

	Dmrg::MatrixVectorStored<ComplexType, FakeTypes> matrixVector(model, connection, aux);

	CHECK(matrixVector.rows() == 2);
	CHECK(matrixVector.cols() == 2);
}

TEST_CASE("MatrixVectorStored applies the stored matrix", "[MatrixVectorStored]")
{
	FakeModel                                        model;
	FakeHamiltonianConnection                        connection(makeMatrix());
	FakeAux                                          aux;
	Dmrg::MatrixVectorStored<ComplexType, FakeTypes> matrixVector(model, connection, aux);

	const std::vector<ComplexType> y { ComplexType(2.0, -1.0), ComplexType(-1.0, 2.0) };
	std::vector<ComplexType>       x(2, ComplexType(0.0, 0.0));
	matrixVector.matrixVectorProduct(x, y);

	checkVector(x, { ComplexType(1.0, 5.0), ComplexType(4.0, 4.0) });
}

TEST_CASE("MatrixVectorStored applies the adjoint product with LdaggerL", "[MatrixVectorStored]")
{
	FakeModel model;
	model.params_.options.ldaggerL = true;
	FakeHamiltonianConnection                        connection(makeMatrix());
	FakeAux                                          aux;
	Dmrg::MatrixVectorStored<ComplexType, FakeTypes> matrixVector(model, connection, aux);

	const std::vector<ComplexType> y { ComplexType(2.0, -1.0), ComplexType(-1.0, 2.0) };
	std::vector<ComplexType>       x(2, ComplexType(0.0, 0.0));
	matrixVector.matrixVectorProduct(x, y);

	checkVector(x, { ComplexType(18.0, -8.0), ComplexType(-2.0, 6.0) });
	CHECK(std::abs(x[0] - ComplexType(-16.0, 18.0)) > 1e-12);
}

TEST_CASE("MatrixVectorStored keeps enabled LdaggerL after model mutation", "[MatrixVectorStored]")
{
	FakeModel model;
	model.params_.options.ldaggerL = true;
	FakeHamiltonianConnection                        connection(makeMatrix());
	FakeAux                                          aux;
	Dmrg::MatrixVectorStored<ComplexType, FakeTypes> matrixVector(model, connection, aux);

	model.params_.options.ldaggerL = false;
	const std::vector<ComplexType> y { ComplexType(2.0, -1.0), ComplexType(-1.0, 2.0) };
	std::vector<ComplexType>       x(2, ComplexType(0.0, 0.0));
	matrixVector.matrixVectorProduct(x, y);

	checkVector(x, { ComplexType(18.0, -8.0), ComplexType(-2.0, 6.0) });
}

TEST_CASE("MatrixVectorStored keeps disabled LdaggerL after model mutation", "[MatrixVectorStored]")
{
	FakeModel                                        model;
	FakeHamiltonianConnection                        connection(makeMatrix());
	FakeAux                                          aux;
	Dmrg::MatrixVectorStored<ComplexType, FakeTypes> matrixVector(model, connection, aux);

	model.params_.options.ldaggerL = true;
	const std::vector<ComplexType> y { ComplexType(2.0, -1.0), ComplexType(-1.0, 2.0) };
	std::vector<ComplexType>       x(2, ComplexType(0.0, 0.0));
	matrixVector.matrixVectorProduct(x, y);

	checkVector(x, { ComplexType(1.0, 5.0), ComplexType(4.0, 4.0) });
}
