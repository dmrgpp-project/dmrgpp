#include "OneTimeEvolution.hpp"
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

TEST_CASE("OneTimeEvolution parses and owns its creation parameters", "[OneTimeEvolution]")
{
	using EvolutionType = Dmrg::OneTimeEvolution<FakePvectors>;

	FakePvectors                       pvectors;
	FakePvectors::VectorWithOffsetType source = { 1.0 };
	EvolutionType                      evolution(
            "|P0>",
            "TimeEvolve{steps=5,tau=0.1,advanceEach=4,disposition=2,algorithm=Chebyshev:1:2,depends=P1:P3,EnergyForExp=2.5}",
            8,
            1.5,
            0.5);
	evolution.initialize(4, source, pvectors);

	CHECK(evolution.indices() == std::vector<SizeType> { 4, 10, 11, 12, 13 });
	CHECK(pvectors.createdTimes == std::vector<double>(4, 0.5));
	CHECK(evolution.timeSteps() == 5);
	CHECK(evolution.disposition() == 2);
	CHECK(evolution.tau() == 0.1);
	CHECK(evolution.advanceEach() == 4);
	CHECK(evolution.algorithm() == "Chebyshev:1:2");
	CHECK(evolution.chebyTransform() == std::vector<double> { 1.0, 2.0 });
	CHECK(evolution.depends() == std::vector<SizeType> { 1, 3 });
	CHECK(evolution.energyForExp() == 2.5);
	CHECK(evolution.sourceToDestroy() == 0);
}

TEST_CASE("OneTimeEvolution does not destroy a non-P source", "[OneTimeEvolution]")
{
	using EvolutionType = Dmrg::OneTimeEvolution<FakePvectors>;

	FakePvectors                       pvectors;
	FakePvectors::VectorWithOffsetType source = { 1.0 };
	EvolutionType evolution("|gs>", "TimeEvolve{steps=2,disposition=2}", 8, 1.5, 0.5);
	evolution.initialize(4, source, pvectors);

	CHECK(evolution.sourceToDestroy() == -1);
}
