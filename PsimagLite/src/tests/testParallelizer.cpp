#include "MatrixMatrix.hpp"
#include <Kokkos_Core.hpp>
#include <PsimagLite/CodeSectionParams.h>
#include <PsimagLite/Parallelizer.h>
#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>

#include <Kokkos_Profiling_ScopedRegion.hpp>

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
Kokkos::Profiling::ScopedRegion region("single");
	MatrixMatrix                           work;
	PsimagLite::Parallelizer<MatrixMatrix> parallelizer(PsimagLite::CodeSectionParams(1));

	REQUIRE(work.tasks() == 10);
  for (int i=0; i<1000; ++i)
    parallelizer.loopCreate(work);
  Kokkos::Timer timer;  
  for (int i=0; i<1000; ++i)                
	  parallelizer.loopCreate(work);
  std::cout << "single: " << timer.seconds() << std::endl;
}



#ifdef USE_PTHREADS
TEST_CASE("Parallelizer concurrently completes independent GEMM tasks", "[Parallelizer][pthread]")
{
Kokkos::Profiling::ScopedRegion region("threaded");

	MatrixMatrix                           work;
	PsimagLite::Parallelizer<MatrixMatrix> parallelizer(PsimagLite::CodeSectionParams(4));

  REQUIRE(work.tasks() == 10);
  for (int i=0; i<1000; ++i)
    parallelizer.loopCreate(work);
  Kokkos::Timer timer;
  for (int i=0; i<1000; ++i)
	  parallelizer.loopCreate(work);
  std::cout << "threaded: " << timer.seconds() << std::endl;
}
#endif

int main(int argc, char* argv[])
{
	Kokkos::ScopeGuard kokkosScopeGuard(argc, argv);
	Catch::Session     session;
	return session.run(argc, argv);
}
