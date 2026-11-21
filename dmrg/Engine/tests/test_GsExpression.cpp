#include "GsExpression.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>

#include <string>
#include <vector>

TEST_CASE("GsExpression parses a PureState expression", "[GsExpression]")
{
	std::vector<SizeType> indices = { 9, 9 };

	Dmrg::GsExpression::setPureState(indices, "PureState[0,2,3,1]");

	CHECK(indices == std::vector<SizeType> { 0, 2, 3, 1 });
}

TEST_CASE("GsExpression leaves indices unchanged for an absent expression", "[GsExpression]")
{
	std::vector<SizeType> indices = { 4, 2 };

	Dmrg::GsExpression::setPureState(indices, "");

	CHECK(indices == std::vector<SizeType> { 4, 2 });
}

TEST_CASE("GsExpression rejects malformed PureState expressions atomically", "[GsExpression]")
{
	std::vector<SizeType> indices = { 7, 8 };

	SECTION("wrong expression form")
	{
		CHECK_THROWS_WITH(Dmrg::GsExpression::setPureState(indices, "PureState[0,1"),
		                  "Expected PureState[...] for gs\n");
	}

	SECTION("non-integer index")
	{
		CHECK_THROWS_WITH(Dmrg::GsExpression::setPureState(indices, "PureState[0,two]"),
		                  "PureState indices must be non-negative integers\n");
	}

	SECTION("negative index")
	{
		CHECK_THROWS_WITH(Dmrg::GsExpression::setPureState(indices, "PureState[-1,0]"),
		                  "PureState indices must be non-negative integers\n");
	}

	CHECK(indices == std::vector<SizeType> { 7, 8 });
}
