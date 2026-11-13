#include <PsimagLite/CanonicalExpression.h>

#include <catch2/catch_test_macros.hpp>

namespace {

struct RecordingAssignAndDestroy {

	using ResultType = PsimagLite::String;

	RecordingAssignAndDestroy(const ResultType& result)
	    : result_(result)
	{ }

	void assignBackward(ResultType& result) const { result = result_; }

	void assign(const RecordingAssignAndDestroy& other) { result_ = other.result_; }

	void plusBackward(ResultType& result) const { result += "+" + result_; }

	void multiply(const RecordingAssignAndDestroy& other) { result_ += "*" + other.result_; }

	void multiplyScalar(const double&) { }

	bool isEmpty() const { return result_.empty(); }

	bool metaEqual(const ResultType&) const { return true; }

	static bool metaEqual(const RecordingAssignAndDestroy&, const RecordingAssignAndDestroy&)
	{
		return true;
	}

private:

	ResultType result_;
};

struct RecordingItemSpec {

	using ResultType        = PsimagLite::String;
	using ComplexOrRealType = double;
	using AuxiliaryType     = PsimagLite::Vector<PsimagLite::String>::Type;

	ResultType operator()(const PsimagLite::String& factor, AuxiliaryType& factors) const
	{
		factors.push_back(factor);
		return factor;
	}

	static bool isEmpty(const ResultType& result) { return result.empty(); }

	static bool metaEqual(const ResultType&, const ResultType&) { return true; }
};

using CanonicalExpressionType
    = PsimagLite::CanonicalExpression<RecordingItemSpec, RecordingAssignAndDestroy>;

} // namespace

TEST_CASE("CanonicalExpression keeps braced operator arguments intact", "[CanonicalExpression]")
{
	RecordingItemSpec::AuxiliaryType    factors;
	RecordingItemSpec                   itemSpec;
	CanonicalExpressionType             canonicalExpression(itemSpec);
	RecordingItemSpec::ResultType       result;
	const RecordingItemSpec::ResultType empty;

	const PsimagLite::String expression = "TimeEvolve{tau=-0.1,formula=a*b+c}*|P0>+|P1>";

	REQUIRE(canonicalExpression(result, expression, empty, factors));
	REQUIRE(factors.size() == 3);
	CHECK(factors[0] == "TimeEvolve{tau=-0.1,formula=a*b+c}");
	CHECK(factors[1] == "|P0>");
	CHECK(factors[2] == "|P1>");
	CHECK(result == "TimeEvolve{tau=-0.1,formula=a*b+c}*|P0>+|P1>");
}

TEST_CASE("CanonicalExpression rejects unmatched braces", "[CanonicalExpression]")
{
	RecordingItemSpec::AuxiliaryType    factors;
	RecordingItemSpec                   itemSpec;
	CanonicalExpressionType             canonicalExpression(itemSpec);
	RecordingItemSpec::ResultType       result;
	const RecordingItemSpec::ResultType empty;

	CHECK_THROWS_AS(canonicalExpression(result, "TimeEvolve{tau=-0.1*|P0>", empty, factors),
	                PsimagLite::RuntimeError);
	CHECK_THROWS_AS(canonicalExpression(result, "TimeEvolve}*|P0>", empty, factors),
	                PsimagLite::RuntimeError);
}

TEST_CASE("CanonicalExpression preserves ordinary factor and term splitting",
          "[CanonicalExpression]")
{
	RecordingItemSpec::AuxiliaryType    factors;
	RecordingItemSpec                   itemSpec;
	CanonicalExpressionType             canonicalExpression(itemSpec);
	RecordingItemSpec::ResultType       result;
	const RecordingItemSpec::ResultType empty;

	REQUIRE(canonicalExpression(result, "A*|P0>+B*|P1>", empty, factors));
	REQUIRE(factors.size() == 4);
	CHECK(factors[0] == "A");
	CHECK(factors[1] == "|P0>");
	CHECK(factors[2] == "B");
	CHECK(factors[3] == "|P1>");
}
