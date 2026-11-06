#include <PsimagLite/Io/IoSelector.h>
#include <PsimagLite/PsimagLite.h>

#include "ProgramGlobals.h"
#include "Wft/WftOptions.h"
#include "Wft/WftParallelismPolicy.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <vector>

namespace {

struct Options {

	bool isSet(const PsimagLite::String&) const { return false; }
};

using WftOptionsType = Dmrg::WftOptions<std::vector<double>, Options>;
using PolicyType     = Dmrg::WftParallelismPolicy;

PsimagLite::CodeSectionParams params(SizeType threads, const PsimagLite::String& name)
{
	PsimagLite::CodeSectionParams result(threads);
	result.name = name;
	return result;
}

} // namespace

TEST_CASE("WFT blocks rejects concurrent outer and inner workers", "[WftParallelismPolicy]")
{
	const auto outer            = params(2, "ParallelWftMany");
	const auto blocks           = params(2, "WftAccelBlocks");
	const auto patchesCtor      = params(1, "BlockDiagWfCtor");
	const auto patchesTransform = params(1, "BlockDiagWfTransform");

	CHECK_THROWS_WITH(
	    PolicyType::check<WftOptionsType>(
	        WftOptionsType::ACCEL_BLOCKS, outer, blocks, patchesCtor, patchesTransform),
	    "ParallelWftMany=2 conflicts with accelerated WFT blocks section "
	    "WftAccelBlocks=2; set one of these thread counts to 1\n");
}

TEST_CASE("WFT patches rejects concurrent outer and constructor workers", "[WftParallelismPolicy]")
{
	const auto outer            = params(2, "ParallelWftMany");
	const auto blocks           = params(1, "WftAccelBlocks");
	const auto patchesCtor      = params(2, "BlockDiagWfCtor");
	const auto patchesTransform = params(1, "BlockDiagWfTransform");

	CHECK_THROWS_WITH(
	    PolicyType::check<WftOptionsType>(
	        WftOptionsType::ACCEL_PATCHES, outer, blocks, patchesCtor, patchesTransform),
	    "ParallelWftMany=2 conflicts with accelerated WFT patches section "
	    "BlockDiagWfCtor=2; set one of these thread counts to 1\n");
}

TEST_CASE("WFT patches rejects concurrent outer and transform workers", "[WftParallelismPolicy]")
{
	const auto outer            = params(2, "ParallelWftMany");
	const auto blocks           = params(1, "WftAccelBlocks");
	const auto patchesCtor      = params(1, "BlockDiagWfCtor");
	const auto patchesTransform = params(2, "BlockDiagWfTransform");

	CHECK_THROWS_WITH(
	    PolicyType::check<WftOptionsType>(
	        WftOptionsType::ACCEL_PATCHES, outer, blocks, patchesCtor, patchesTransform),
	    "ParallelWftMany=2 conflicts with accelerated WFT patches section "
	    "BlockDiagWfTransform=2; set one of these thread counts to 1\n");
}

TEST_CASE("WFT parallelism policy permits non-conflicting configurations", "[WftParallelismPolicy]")
{
	const auto outerOne     = params(1, "ParallelWftMany");
	const auto outerTwo     = params(2, "ParallelWftMany");
	const auto oneBlock     = params(1, "WftAccelBlocks");
	const auto twoBlocks    = params(2, "WftAccelBlocks");
	const auto oneCtor      = params(1, "BlockDiagWfCtor");
	const auto twoCtor      = params(2, "BlockDiagWfCtor");
	const auto oneTransform = params(1, "BlockDiagWfTransform");
	const auto twoTransform = params(2, "BlockDiagWfTransform");

	CHECK_NOTHROW(PolicyType::check<WftOptionsType>(
	    WftOptionsType::ACCEL_BLOCKS, outerOne, twoBlocks, oneCtor, oneTransform));
	CHECK_NOTHROW(PolicyType::check<WftOptionsType>(
	    WftOptionsType::ACCEL_BLOCKS, outerTwo, oneBlock, oneCtor, oneTransform));
	CHECK_NOTHROW(PolicyType::check<WftOptionsType>(
	    WftOptionsType::ACCEL_PATCHES, outerTwo, oneBlock, oneCtor, oneTransform));
	CHECK_NOTHROW(PolicyType::check<WftOptionsType>(
	    WftOptionsType::ACCEL_NONE, outerTwo, twoBlocks, twoCtor, twoTransform));
	CHECK_NOTHROW(PolicyType::check<WftOptionsType>(
	    WftOptionsType::ACCEL_SVD, outerTwo, twoBlocks, twoCtor, twoTransform));
}
