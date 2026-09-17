#include "GroupOfOneTimeEvolutions.h"
#include <catch2/catch_test_macros.hpp>
#include <string>
#include <utility>
#include <vector>

namespace {

class FakePvectors {

public:

	using RealType             = double;
	using VectorWithOffsetType = std::vector<double>;

	template <typename SomeLambdaType>
	void createNew(const VectorWithOffsetType&, RealType time, SomeLambdaType& lambda)
	{
		createdTimes.push_back(time);
		lambda(nextIndex++);
	}

	void registerLastKrylovSlot(SizeType pIndex, SizeType slotIndex)
	{
		lastKrylovSlot = { pIndex, slotIndex };
	}

	SizeType                      nextIndex = 10;
	std::vector<RealType>         createdTimes;
	std::pair<SizeType, SizeType> lastKrylovSlot = { 0, 0 };
};

} // namespace

TEST_CASE("Time evolution inherits its nonzero initial time", "[GroupOfOneTimeEvolutions]")
{
	using GroupType     = Dmrg::GroupOfOneTimeEvolutions<FakePvectors>;
	using EvolutionType = typename GroupType::OneTimeEvolutionType;

	FakePvectors                        pvectors;
	FakePvectors::VectorWithOffsetType  source = { 1.0 };
	EvolutionType                       evolution(4, source, "|P0>", 0, 5, 0.5, pvectors);
	const std::vector<SizeType>         expectedIndices        = { 4, 10, 11, 12, 13 };
	const std::pair<SizeType, SizeType> expectedLastKrylovSlot = { 4, 13 };

	CHECK(evolution.time() == 0.5);
	CHECK(evolution.indices() == expectedIndices);
	CHECK(pvectors.createdTimes == std::vector<double>(4, 0.5));
	CHECK(pvectors.lastKrylovSlot == expectedLastKrylovSlot);

	evolution.advanceTime(0.1);
	CHECK(evolution.time() == 0.6);
}
