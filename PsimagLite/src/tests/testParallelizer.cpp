#include "MatrixMatrix.hpp"
#include <Kokkos_Core.hpp>
#include <PsimagLite/CodeSectionParams.h>
#include <PsimagLite/Concurrency.h>
#include <PsimagLite/Parallelizer.h>
#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>

namespace {

MatrixMatrix::OutputMatrix referenceProduct(const MatrixMatrix::LeftMatrix&  left,
                                            const MatrixMatrix::RightMatrix& right)
{
	MatrixMatrix::OutputMatrix expected {};
	for (SizeType column = 0; column < 3; ++column) {
		for (SizeType row = 0; row < 3; ++row) {
			for (SizeType inner = 0; inner < 8; ++inner)
				expected[row + column * 3]
				    += left[row + inner * 3] * right[inner + column * 8];
		}
	}

	return expected;
}

void checkOutputs(const MatrixMatrix& work)
{
	for (SizeType task = 0; task < work.tasks(); ++task) {
		const auto  expected = referenceProduct(work.left(task), work.right(task));
		const auto& actual   = work.output(task);
		for (SizeType column = 0; column < 3; ++column) {
			for (SizeType row = 0; row < 3; ++row) {
				INFO("task=" << task << ", row=" << row << ", column=" << column);
				CHECK(actual[row + column * 3] == expected[row + column * 3]);
			}
		}
	}
}

} // namespace

TEST_CASE("Concurrency resolves named code section parameters", "[Parallelizer][CodeSectionParams]")
{
	using ConcurrencyType = PsimagLite::Concurrency;

	ConcurrencyType::clearCodeSectionParamsOverrides();
	ConcurrencyType::setOptions(PsimagLite::CodeSectionParams(1));

	const auto kronDefault      = ConcurrencyType::codeSectionParamsFor("KronMatrix");
	const auto operatorsDefault = ConcurrencyType::codeSectionParamsFor("Operators");
	REQUIRE(kronDefault.name == "KronMatrix");
	REQUIRE(kronDefault.npthreads == 1);
	REQUIRE(operatorsDefault.name == "Operators");
	REQUIRE(operatorsDefault.npthreads == 1);

	REQUIRE_THROWS_AS(ConcurrencyType::setCodeSectionParamsOverride("NotASection", 1),
	                  PsimagLite::RuntimeError);
	REQUIRE_THROWS_AS(ConcurrencyType::setCodeSectionParamsOverride("KronMatrix", 0),
	                  PsimagLite::RuntimeError);

#ifdef USE_PTHREADS
	ConcurrencyType::setOptions(PsimagLite::CodeSectionParams(4));
	ConcurrencyType::setCodeSectionParamsOverride("KronMatrix", 2);
	REQUIRE(ConcurrencyType::codeSectionParamsFor("KronMatrix").npthreads == 2);
	REQUIRE(ConcurrencyType::codeSectionParamsFor("Operators").npthreads == 4);
#endif

	ConcurrencyType::clearCodeSectionParamsOverrides();
	ConcurrencyType::setOptions(PsimagLite::CodeSectionParams(1));
}

TEST_CASE("Parallelizer completes independent GEMM tasks", "[Parallelizer]")
{
	MatrixMatrix                           work;
	PsimagLite::Parallelizer<MatrixMatrix> parallelizer(PsimagLite::CodeSectionParams(1));

	REQUIRE(work.tasks() == 4);
	parallelizer.loopCreate(work);

	checkOutputs(work);
}

#ifdef USE_PTHREADS
TEST_CASE("Parallelizer concurrently completes independent GEMM tasks", "[Parallelizer][pthread]")
{
	MatrixMatrix                           work;
	PsimagLite::Parallelizer<MatrixMatrix> parallelizer(PsimagLite::CodeSectionParams(4));

	parallelizer.loopCreate(work);

	checkOutputs(work);
}
#endif

int main(int argc, char* argv[])
{
	Kokkos::ScopeGuard kokkosScopeGuard(argc, argv);
	Catch::Session     session;
	return session.run(argc, argv);
}
