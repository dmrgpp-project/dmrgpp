#include <PsimagLite/GetBraOrKet.h>
#include <PsimagLite/PsimagLite.h>
#include <catch2/catch_test_macros.hpp>

using G = PsimagLite::GetBraOrKet;

TEST_CASE("GetBraOrKet ket vs bra detection", "[GetBraOrKet]")
{
	CHECK(G("|gs>").isKet() == true);
	CHECK(G("<gs|").isKet() == false);
	CHECK(G("|P3>").isKet() == true);
	CHECK(G("<P3|").isKet() == false);
}

TEST_CASE("GetBraOrKet ground state", "[GetBraOrKet]")
{
	G gs("|gs>");
	CHECK(gs.isPvector() == false);
	CHECK(gs.isRvector() == false);
	CHECK(gs.isLastKrylov() == false);
}

TEST_CASE("GetBraOrKet P-vector index parsing", "[GetBraOrKet]")
{
	SECTION("P0")
	{
		G p0("|P0>");
		CHECK(p0.isPvector() == true);
		CHECK(p0.pIndex() == 0);
		CHECK(p0.isLastKrylov() == false);
	}

	SECTION("P2")
	{
		G p2("|P2>");
		CHECK(p2.isPvector() == true);
		CHECK(p2.pIndex() == 2);
		CHECK(p2.isLastKrylov() == false);
	}

	SECTION("P12 multi-digit")
	{
		G p12("|P12>");
		CHECK(p12.isPvector() == true);
		CHECK(p12.pIndex() == 12);
		CHECK(p12.isLastKrylov() == false);
	}

	SECTION("bra form <P2|")
	{
		G p2bra("<P2|");
		CHECK(p2bra.isPvector() == true);
		CHECK(p2bra.pIndex() == 2);
		CHECK(p2bra.isKet() == false);
	}
}

TEST_CASE("GetBraOrKet .last suffix", "[GetBraOrKet]")
{
	SECTION("ket |P2.last>")
	{
		G p2last("|P2.last>");
		CHECK(p2last.isPvector() == true);
		CHECK(p2last.pIndex() == 2);
		CHECK(p2last.isLastKrylov() == true);
		CHECK(p2last.isKet() == true);
	}

	SECTION("bra <P2.last|")
	{
		G p2last_bra("<P2.last|");
		CHECK(p2last_bra.isPvector() == true);
		CHECK(p2last_bra.pIndex() == 2);
		CHECK(p2last_bra.isLastKrylov() == true);
		CHECK(p2last_bra.isKet() == false);
	}

	SECTION("P0.last")
	{
		G p0last("|P0.last>");
		CHECK(p0last.pIndex() == 0);
		CHECK(p0last.isLastKrylov() == true);
	}

	SECTION("regular P2 does not have .last set")
	{
		G p2("|P2>");
		CHECK(p2.isLastKrylov() == false);
	}
}
