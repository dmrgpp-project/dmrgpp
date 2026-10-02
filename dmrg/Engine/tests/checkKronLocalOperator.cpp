#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using ComplexType = std::complex<double>;
using MatrixType  = std::vector<std::vector<ComplexType>>;

bool isFinite(const ComplexType& value)
{
	return std::isfinite(value.real()) && std::isfinite(value.imag());
}

std::vector<std::string> readLines(const std::string& filename)
{
	std::ifstream input(filename);
	if (!input)
		throw std::runtime_error("Cannot open " + filename);

	std::vector<std::string> lines;
	std::string              line;
	while (std::getline(input, line))
		lines.push_back(line);
	return lines;
}

std::string valueAfter(const std::string& line, const std::string& prefix)
{
	if (line.compare(0, prefix.size(), prefix) != 0)
		throw std::runtime_error("Expected " + prefix + " line, found: " + line);
	return line.substr(prefix.size());
}

ComplexType parseComplex(const std::string& text)
{
	std::istringstream input(text);
	ComplexType        value;
	input >> value;
	if (!input || !isFinite(value))
		throw std::runtime_error("Invalid complex value: " + text);
	input >> std::ws;
	if (input.peek() != std::char_traits<char>::eof())
		throw std::runtime_error("Invalid complex value: " + text);
	return value;
}

std::size_t parseSize(const std::string& text, const std::string& description)
{
	std::istringstream input(text);
	std::size_t        value = 0;
	input >> value;
	if (!input)
		throw std::runtime_error("Invalid " + description + ": " + text);
	input >> std::ws;
	if (input.peek() != std::char_traits<char>::eof())
		throw std::runtime_error("Invalid " + description + ": " + text);
	return value;
}

bool columnsClose(const std::vector<ComplexType>& direct,
                  const MatrixType&               matrix,
                  std::size_t                     column,
                  double                          tolerance)
{
	if (direct.size() != matrix.size() || column >= matrix.size())
		return false;
	for (std::size_t row = 0; row < direct.size(); ++row) {
		if (!isFinite(direct[row]) || !isFinite(matrix[row][column])
		    || std::abs(direct[row] - matrix[row][column]) > tolerance)
			return false;
	}
	return true;
}

void verifyOrientationCheck()
{
	const MatrixType               matrix       = { { 0, 1 }, { 2, 0 } };
	const std::vector<ComplexType> directColumn = { 0, 2 };
	if (!columnsClose(directColumn, matrix, 0, 0))
		throw std::runtime_error("Orientation-check self-test rejected a correct column");
	MatrixType transposed = { { 0, 2 }, { 1, 0 } };
	if (columnsClose(directColumn, transposed, 0, 0))
		throw std::runtime_error("Orientation-check self-test accepted a transpose");
}

void requireClose(const ComplexType& actual,
                  const ComplexType& expected,
                  double             tolerance,
                  const std::string& message)
{
	if (!isFinite(actual) || !isFinite(expected) || std::abs(actual - expected) > tolerance) {
		std::ostringstream output;
		output << message << ": actual=" << actual << " expected=" << expected;
		throw std::runtime_error(output.str());
	}
}

void verify(const std::string& filename)
{
	const std::vector<std::string> lines = readLines(filename);
	if (std::count(lines.begin(), lines.end(), "LocalOperator") != 1)
		throw std::runtime_error("Expected exactly one LocalOperator section");

	auto line = std::find(lines.begin(), lines.end(), "LocalOperator");
	++line;
	if (std::distance(line, lines.end()) < 10)
		throw std::runtime_error("LocalOperator section is truncated");
	if (valueAfter(*line++, "OperatorName=") != "sz")
		throw std::runtime_error("Expected the sz operator");
	if (valueAfter(*line++, "OperatorSite=") != "5")
		throw std::runtime_error("Expected operator site 5");
	if (valueAfter(*line++, "Rows=") != "70 Cols=70")
		throw std::runtime_error("Expected a 70 x 70 operator");
	if (valueAfter(*line++, "BasisOrder=") != "KronPatchPacked")
		throw std::runtime_error("Unexpected operator basis order");
	if (valueAfter(*line++, "Direction=") != "EXPAND_ENVIRON")
		throw std::runtime_error("Unexpected sweep direction");
	if (valueAfter(*line++, "Border=") != "0")
		throw std::runtime_error("Unexpected border mode");
	if (valueAfter(*line++, "MatrixAction=") != "dest=O*src")
		throw std::runtime_error("Unexpected matrix action");
	if (valueAfter(*line++, "MatrixFormat=") != "matrix-market-coordinate")
		throw std::runtime_error("Unexpected matrix format");
	if (*line++ != "%%MatrixMarket matrix coordinate real general")
		throw std::runtime_error("Expected a real Matrix Market payload");

	std::size_t rows = 0;
	std::size_t cols = 0;
	std::size_t nnz  = 0;
	{
		std::istringstream dimensions(*line++);
		dimensions >> rows >> cols >> nnz;
		if (!dimensions || rows != 70 || cols != 70)
			throw std::runtime_error("Invalid Matrix Market dimensions");
		dimensions >> std::ws;
		if (dimensions.peek() != std::char_traits<char>::eof())
			throw std::runtime_error("Invalid Matrix Market dimensions");
	}

	if (nnz > static_cast<std::size_t>(std::distance(line, lines.end())))
		throw std::runtime_error("Matrix Market payload is truncated");
	MatrixType matrix(rows, std::vector<ComplexType>(cols, 0));
	for (std::size_t entry = 0; entry < nnz; ++entry, ++line) {
		std::size_t        row   = 0;
		std::size_t        col   = 0;
		double             value = 0;
		std::istringstream input(*line);
		input >> row >> col >> value;
		if (!input || !std::isfinite(value) || row == 0 || row > rows || col == 0
		    || col > cols)
			throw std::runtime_error("Invalid Matrix Market entry: " + *line);
		input >> std::ws;
		if (input.peek() != std::char_traits<char>::eof())
			throw std::runtime_error("Invalid Matrix Market entry: " + *line);
		matrix[row - 1][col - 1] += value;
	}

	const double tolerance = 1e-12;
	for (std::size_t row = 0; row < rows; ++row) {
		for (std::size_t col = 0; col < cols; ++col) {
			requireClose(matrix[row][col],
			             std::conj(matrix[col][row]),
			             tolerance,
			             "Local operator is not Hermitian");
		}
	}

	line = std::find(line, lines.end(), "Vector");
	if (line == lines.end())
		throw std::runtime_error("Vector section is missing");
	++line;
	if (line == lines.end())
		throw std::runtime_error("Vector size is missing");
	const std::size_t vectorSize = parseSize(*line++, "vector size");
	if (vectorSize != rows)
		throw std::runtime_error("Operator and vector dimensions differ");
	if (vectorSize > static_cast<std::size_t>(std::distance(line, lines.end())))
		throw std::runtime_error("Vector section is truncated");

	std::vector<ComplexType> vector(vectorSize);
	for (ComplexType& value : vector)
		value = parseComplex(*line++);

	const std::size_t sampleColumns[] = { 0, cols / 2, cols - 1 };
	for (const std::size_t column : sampleColumns) {
		if (line == lines.end()
		    || valueAfter(*line++, "OperatorDirectColumn=") != std::to_string(column))
			throw std::runtime_error("Expected direct operator column "
			                         + std::to_string(column));
		if (line == lines.end())
			throw std::runtime_error("Direct operator column size is missing");
		const std::size_t directSize = parseSize(*line++, "direct operator column size");
		if (directSize != rows
		    || directSize > static_cast<std::size_t>(std::distance(line, lines.end())))
			throw std::runtime_error("Direct operator column has an invalid size");
		std::vector<ComplexType> directColumn(directSize);
		for (ComplexType& value : directColumn)
			value = parseComplex(*line++);
		if (!columnsClose(directColumn, matrix, column, tolerance))
			throw std::runtime_error(
			    "Serialized matrix column disagrees with direct application");
	}

	std::vector<ComplexType> result(rows, 0);
	for (std::size_t row = 0; row < rows; ++row)
		for (std::size_t col = 0; col < cols; ++col)
			result[row] += matrix[row][col] * vector[col];

	ComplexType serializedExpectation = 0;
	for (std::size_t i = 0; i < vectorSize; ++i)
		serializedExpectation += std::conj(vector[i]) * result[i];

	line = std::find_if(line,
	                    lines.end(),
	                    [](const std::string& text)
	                    { return text.compare(0, 20, "OperatorExpectation=") == 0; });
	if (line == lines.end())
		throw std::runtime_error("OperatorExpectation is missing");
	std::istringstream expectation(valueAfter(*line, "OperatorExpectation="));
	double             realPart = 0;
	double             imagPart = 0;
	expectation >> realPart >> imagPart;
	if (!expectation || !std::isfinite(realPart) || !std::isfinite(imagPart))
		throw std::runtime_error("Invalid OperatorExpectation");
	expectation >> std::ws;
	if (expectation.peek() != std::char_traits<char>::eof())
		throw std::runtime_error("Invalid OperatorExpectation");
	requireClose(serializedExpectation,
	             ComplexType(realPart, imagPart),
	             tolerance,
	             "Serialized and direct-apply expectations differ");
}

} // namespace

int main(int argc, char* argv[])
{
	if (argc != 2) {
		std::cerr << "Usage: " << argv[0] << " KRON_DUMP\n";
		return EXIT_FAILURE;
	}

	try {
		verifyOrientationCheck();
		verify(argv[1]);
	} catch (const std::exception& error) {
		std::cerr << error.what() << '\n';
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
