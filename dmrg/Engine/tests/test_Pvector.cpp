#include <PsimagLite/PsimagLite.h>

#include "Pvector.h"
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <complex>
#include <limits>

TEMPLATE_TEST_CASE("Pvector stores physical time", "[Pvector]", double, std::complex<double>)
{
	Dmrg::Pvector<TestType> pvector("|P0>");
	CHECK(pvector.time() == 0.0);

	pvector.setTime(0.25);
	const Dmrg::Pvector<TestType>& constPvector = pvector;
	CHECK(constPvector.time() == 0.25);
	CHECK(constPvector.weight() == 1.0);
	CHECK(constPvector.lastName() == "|P0>");
	CHECK_FALSE(constPvector.isDone());
}

TEMPLATE_TEST_CASE("Pvector sums retain equal physical time",
                   "[Pvector]",
                   double,
                   std::complex<double>)
{
	Dmrg::Pvector<TestType> lhs("|P0>");
	Dmrg::Pvector<TestType> rhs("|P1>");
	lhs.setTime(0.25);
	rhs.setTime(0.25);
	lhs.setAsDone();

	lhs.sum(rhs, "|P0>+|P1>");

	CHECK(lhs.isDone());
	CHECK(lhs.time() == 0.25);
}

TEMPLATE_TEST_CASE("Pvector sums times equal within floating-point roundoff",
                   "[Pvector]",
                   double,
                   std::complex<double>)
{
	using RealType = typename PsimagLite::Real<TestType>::Type;

	Dmrg::Pvector<TestType> lhs("|P0>");
	Dmrg::Pvector<TestType> rhs("|P1>");
	lhs.setTime(RealType(0.1) + RealType(0.2));
	rhs.setTime(RealType(0.3));
	lhs.setAsDone();

	CHECK(lhs.hasSameTime(rhs));
	CHECK_NOTHROW(lhs.sum(rhs, "|P0>+|P1>"));
}

TEMPLATE_TEST_CASE("Pvector rejects times beyond floating-point roundoff",
                   "[Pvector]",
                   double,
                   std::complex<double>)
{
	using RealType = typename PsimagLite::Real<TestType>::Type;

	Dmrg::Pvector<TestType> lhs("|P0>");
	Dmrg::Pvector<TestType> rhs("|P1>");
	lhs.setTime(RealType(0.25));
	rhs.setTime(RealType(0.25) + RealType(1e-8));
	lhs.setAsDone();

	CHECK_FALSE(lhs.hasSameTime(rhs));
	CHECK_THROWS(lhs.sum(rhs, "|P0>+|P1>"));
}

TEMPLATE_TEST_CASE("Pvector sums equal infinite times", "[Pvector]", double, std::complex<double>)
{
	using RealType = typename PsimagLite::Real<TestType>::Type;

	Dmrg::Pvector<TestType> lhs("|P0>");
	Dmrg::Pvector<TestType> rhs("|P1>");
	lhs.setTime(std::numeric_limits<RealType>::infinity());
	rhs.setTime(std::numeric_limits<RealType>::infinity());
	lhs.setAsDone();

	CHECK(lhs.hasSameTime(rhs));
	CHECK_NOTHROW(lhs.sum(rhs, "|P0>+|P1>"));
}

TEMPLATE_TEST_CASE("Pvector rejects unequal non-finite times",
                   "[Pvector]",
                   double,
                   std::complex<double>)
{
	using RealType = typename PsimagLite::Real<TestType>::Type;

	Dmrg::Pvector<TestType> lhs("|P0>");
	Dmrg::Pvector<TestType> rhs("|P1>");
	lhs.setTime(std::numeric_limits<RealType>::infinity());
	rhs.setTime(RealType(0));
	lhs.setAsDone();

	CHECK_FALSE(lhs.hasSameTime(rhs));
	CHECK_THROWS(lhs.sum(rhs, "|P0>+|P1>"));

	rhs.setTime(-std::numeric_limits<RealType>::infinity());
	CHECK_FALSE(lhs.hasSameTime(rhs));
	CHECK_THROWS(lhs.sum(rhs, "|P0>+|P1>"));
}

TEMPLATE_TEST_CASE("Pvector rejects sums at different physical times",
                   "[Pvector]",
                   double,
                   std::complex<double>)
{
	Dmrg::Pvector<TestType> lhs("|P0>");
	Dmrg::Pvector<TestType> rhs("|P1>");
	lhs.setTime(0.25);
	rhs.setTime(0.5);
	lhs.setAsDone();
	const SizeType sizeBefore = lhs.size();

	CHECK_THROWS(lhs.sum(rhs, "|P0>+|P1>"));
	CHECK(lhs.time() == 0.25);
	CHECK(lhs.size() == sizeBefore);
	CHECK(lhs.firstName() == "|P0>");
	CHECK(lhs.lastName() == "DONE");
}
