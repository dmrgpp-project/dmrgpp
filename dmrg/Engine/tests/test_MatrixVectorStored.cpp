#include "MatrixVectorStored.h"
#include <PsimagLite/CrsMatrix.h>
#include <PsimagLite/Matrix.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
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
	enum class MatrixSolverEnum
	{
		DENSE,
		LANCZOS,
		ARNOLDISAI
	};

	FakeOptions      options;
	int              maxMatrixRankStored = 10;
	MatrixSolverEnum matrix_solver_enum  = MatrixSolverEnum::LANCZOS;
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

using MatrixVectorType = Dmrg::MatrixVectorStored<ComplexType, FakeTypes>;

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

TEST_CASE("MatrixVectorStored accumulates the stored matrix product", "[MatrixVectorStored]")
{
	FakeModel                                        model;
	FakeHamiltonianConnection                        connection(makeMatrix());
	FakeAux                                          aux;
	Dmrg::MatrixVectorStored<ComplexType, FakeTypes> matrixVector(model, connection, aux);

	const std::vector<ComplexType> y { ComplexType(2.0, -1.0), ComplexType(-1.0, 2.0) };
	std::vector<ComplexType>       x { ComplexType(5.0, 2.0), ComplexType(-3.0, 4.0) };
	matrixVector.matrixVectorProduct(x, y);

	checkVector(x, { ComplexType(6.0, 7.0), ComplexType(1.0, 8.0) });
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

TEST_CASE("MatrixVectorStored accumulates the adjoint product with LdaggerL",
          "[MatrixVectorStored]")
{
	FakeModel model;
	model.params_.options.ldaggerL = true;
	FakeHamiltonianConnection                        connection(makeMatrix());
	FakeAux                                          aux;
	Dmrg::MatrixVectorStored<ComplexType, FakeTypes> matrixVector(model, connection, aux);

	const std::vector<ComplexType> y { ComplexType(2.0, -1.0), ComplexType(-1.0, 2.0) };
	std::vector<ComplexType>       x { ComplexType(5.0, 2.0), ComplexType(-3.0, 4.0) };
	matrixVector.matrixVectorProduct(x, y);

	checkVector(x, { ComplexType(23.0, -6.0), ComplexType(-5.0, 10.0) });
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

TEST_CASE("MatrixVectorStored returns its CRS matrix without LdaggerL", "[MatrixVectorStored]")
{
	FakeModel                 model;
	FakeHamiltonianConnection connection(makeMatrix());
	FakeAux                   aux;
	MatrixVectorType          matrixVector(model, connection, aux);

	const std::vector<ComplexType> y { ComplexType(2.0, -1.0), ComplexType(-1.0, 2.0) };
	std::vector<ComplexType>       x(2, ComplexType(0.0, 0.0));
	matrixVector.toCRS().matrixVectorProduct(x, y);

	checkVector(x, { ComplexType(1.0, 5.0), ComplexType(4.0, 4.0) });
}

TEST_CASE("MatrixVectorStored rejects CRS access with LdaggerL", "[MatrixVectorStored]")
{
	FakeModel model;
	model.params_.options.ldaggerL = true;
	FakeHamiltonianConnection connection(makeMatrix());
	FakeAux                   aux;
	MatrixVectorType          matrixVector(model, connection, aux);

	REQUIRE_THROWS_WITH(matrixVector.toCRS(),
	                    "MatrixVectorStored::toCRS is unavailable with LdaggerL\n");
}

TEST_CASE("MatrixVectorStored rejects dense diagonalization with LdaggerL", "[MatrixVectorStored]")
{
	FakeModel model;
	model.params_.options.ldaggerL = true;
	FakeHamiltonianConnection       connection(makeMatrix());
	FakeAux                         aux;
	MatrixVectorType                matrixVector(model, connection, aux);
	std::vector<double>             eigs;
	PsimagLite::Matrix<ComplexType> eigenvectors;

	REQUIRE_THROWS_WITH(matrixVector.fullDiag(eigs, eigenvectors),
	                    "MatrixVectorStored::fullDiag is unavailable with LdaggerL\n");
}

TEST_CASE("MatrixVectorStored supports LdaggerL with Lanczos", "[MatrixVectorStored]")
{
	FakeModel model;
	model.params_.options.ldaggerL   = true;
	model.params_.matrix_solver_enum = FakeParameters::MatrixSolverEnum::LANCZOS;
	FakeHamiltonianConnection connection(makeMatrix());
	FakeAux                   aux;

	REQUIRE_NOTHROW(MatrixVectorType(model, connection, aux));
}

TEST_CASE("MatrixVectorStored rejects LdaggerL with Dense", "[MatrixVectorStored]")
{
	FakeModel model;
	model.params_.options.ldaggerL   = true;
	model.params_.matrix_solver_enum = FakeParameters::MatrixSolverEnum::DENSE;
	FakeHamiltonianConnection connection(makeMatrix());
	FakeAux                   aux;

	REQUIRE_THROWS_WITH(
	    MatrixVectorType(model, connection, aux),
	    "LdaggerL with MatrixVectorStored does not support MatrixSolver=Dense\n");
}

TEST_CASE("MatrixVectorStored rejects LdaggerL with ArnoldiSaI", "[MatrixVectorStored]")
{
	FakeModel model;
	model.params_.options.ldaggerL   = true;
	model.params_.matrix_solver_enum = FakeParameters::MatrixSolverEnum::ARNOLDISAI;
	FakeHamiltonianConnection connection(makeMatrix());
	FakeAux                   aux;

	REQUIRE_THROWS_WITH(
	    MatrixVectorType(model, connection, aux),
	    "LdaggerL with MatrixVectorStored does not support MatrixSolver=ArnoldiSaI\n");
}

TEST_CASE("MatrixVectorStored supports every solver without LdaggerL", "[MatrixVectorStored]")
{
	const FakeParameters::MatrixSolverEnum solvers[]
	    = { FakeParameters::MatrixSolverEnum::DENSE,
		FakeParameters::MatrixSolverEnum::LANCZOS,
		FakeParameters::MatrixSolverEnum::ARNOLDISAI };

	for (const auto solver : solvers) {
		FakeModel model;
		model.params_.matrix_solver_enum = solver;
		FakeHamiltonianConnection connection(makeMatrix());
		FakeAux                   aux;

		REQUIRE_NOTHROW(MatrixVectorType(model, connection, aux));
	}
}
