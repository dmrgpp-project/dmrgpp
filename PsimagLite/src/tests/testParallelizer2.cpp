#include <Kokkos_Core.hpp>
#include <PsimagLite/BLAS.h>
#include <PsimagLite/CodeSectionParams.h>
#include <PsimagLite/Parallelizer2.h>
#include <array>
#include <atomic>
#include <cassert>
#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>

namespace {

constexpr SizeType FirstTask = 5;
constexpr SizeType Tasks     = 4;

struct MatrixData {
	using LeftMatrix   = std::array<double, 3 * 8>;
	using RightMatrix  = std::array<double, 8 * 3>;
	using OutputMatrix = std::array<double, 3 * 3>;

	MatrixData()
	{
		for (SizeType task = 0; task < Tasks; ++task) {
			for (SizeType column = 0; column < 8; ++column) {
				for (SizeType row = 0; row < 3; ++row)
					left[task][row + column * 3] = static_cast<double>(
					    100 * task + 10 * column + row + 1);
			}

			for (SizeType column = 0; column < 3; ++column) {
				for (SizeType row = 0; row < 8; ++row)
					right[task][row + column * 8]
					    = static_cast<double>(50 * task + 8 * column + row + 1);
			}

			output[task].fill(-999.0);
			calls[task].store(0);
		}
	}

	std::array<LeftMatrix, Tasks>            left;
	std::array<RightMatrix, Tasks>           right;
	std::array<OutputMatrix, Tasks>          output;
	std::array<std::atomic<SizeType>, Tasks> calls;
};

MatrixData::OutputMatrix referenceProduct(const MatrixData::LeftMatrix&  left,
                                          const MatrixData::RightMatrix& right)
{
	MatrixData::OutputMatrix expected {};
	for (SizeType column = 0; column < 3; ++column) {
		for (SizeType row = 0; row < 3; ++row) {
			for (SizeType inner = 0; inner < 8; ++inner)
				expected[row + column * 3]
				    += left[row + inner * 3] * right[inner + column * 8];
		}
	}

	return expected;
}

auto lambdaFor(MatrixData& work)
{
	return [&work](SizeType taskNumber, SizeType)
	{
		assert(taskNumber >= FirstTask);
		assert(taskNumber < FirstTask + Tasks);
		const SizeType task = taskNumber - FirstTask;
		work.calls[task].fetch_add(1, std::memory_order_relaxed);
		psimag::BLAS::GEMM('N',
		                   'N',
		                   3,
		                   3,
		                   8,
		                   1.0,
		                   work.left[task].data(),
		                   3,
		                   work.right[task].data(),
		                   8,
		                   0.0,
		                   work.output[task].data(),
		                   3);
	};
}

void checkOutputs(const MatrixData& work)
{
	for (SizeType task = 0; task < Tasks; ++task) {
		const auto  expected = referenceProduct(work.left[task], work.right[task]);
		const auto& actual   = work.output[task];
		REQUIRE(work.calls[task].load() == 1);
		for (SizeType column = 0; column < 3; ++column) {
			for (SizeType row = 0; row < 3; ++row) {
				INFO("task=" << task << ", row=" << row << ", column=" << column);
				CHECK(actual[row + column * 3] == expected[row + column * 3]);
			}
		}
	}
}

} // namespace

TEST_CASE("Parallelizer2 completes independent GEMM tasks", "[Parallelizer2]")
{
	MatrixData                  work;
	PsimagLite::Parallelizer2<> parallelizer(PsimagLite::CodeSectionParams(1));

	parallelizer.parallelFor(FirstTask, FirstTask + Tasks, lambdaFor(work));

	checkOutputs(work);
}

#if defined(USE_PTHREADS) || defined(_OPENMP)
TEST_CASE("Parallelizer2 concurrently completes independent GEMM tasks",
          "[Parallelizer2][threaded]")
{
	MatrixData                  work;
	PsimagLite::Parallelizer2<> parallelizer(PsimagLite::CodeSectionParams(4));

	parallelizer.parallelFor(FirstTask, FirstTask + Tasks, lambdaFor(work));

	checkOutputs(work);
}
#endif

int main(int argc, char* argv[])
{
	Kokkos::ScopeGuard kokkosScopeGuard(argc, argv);
	Catch::Session     session;
	return session.run(argc, argv);
}
