#include "MatrixMatrix.hpp"
#include <Kokkos_Core.hpp>
#include <PsimagLite/CodeSectionParams.h>
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
