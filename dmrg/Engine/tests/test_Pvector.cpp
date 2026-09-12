#include <PsimagLite/PsimagLite.h>

#include "Pvector.h"
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <complex>

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
