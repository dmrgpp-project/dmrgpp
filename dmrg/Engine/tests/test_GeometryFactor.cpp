#include "InputCheck.h"
#include "ProgramGlobals.h"
#include <PsimagLite/InputNg.h>
#include <catch2/catch_test_macros.hpp>
#include <complex>

namespace {
using InputNgType = PsimagLite::InputNg<Dmrg::InputCheck>;
using GeometryType
    = PsimagLite::Geometry<std::complex<double>, InputNgType::Readable, Dmrg::ProgramGlobals>;

PsimagLite::String twoTermInput()
{
	return R"(##Ainur1.0
TotalNumberOfSites=8;
NumberOfTerms=2;
gt0:GeometryKind="ladder";
gt0:GeometryOptions="ConstantValues";
gt0:LadderLeg=2;
gt0:dir0:Connectors=[1.0];
gt0:dir1:Connectors=[1.0];
gt0:dir0:GeometryFactor=11;
gt0:dir1:GeometryFactor=12;
gt1:GeometryKind="ladder";
gt1:GeometryOptions="ConstantValues";
gt1:LadderLeg=2;
gt1:dir0:Connectors=[1.0];
gt1:dir1:Connectors=[1.0];
gt1:dir0:GeometryFactor=21;
gt1:dir1:GeometryFactor=22;
)";
}

PsimagLite::String oneTermInput(const PsimagLite::String& factor)
{
	return "##Ainur1.0\n"
	       "TotalNumberOfSites=8;\n"
	       "NumberOfTerms=1;\n"
	       "GeometryKind=\"ladder\";\n"
	       "GeometryOptions=\"ConstantValues\";\n"
	       "LadderLeg=2;\n"
	       "dir0:Connectors=[1.0];\n"
	       "dir1:Connectors=[1.0];\n"
	    + factor;
}
} // namespace

TEST_CASE("GeometryFactor is selected by term and direction", "[GeometryFactor]")
{
	Dmrg::InputCheck       inputCheck;
	InputNgType::Writeable writeable(inputCheck, twoTermInput());
	InputNgType::Readable  readable(writeable);
	GeometryType           geometry(readable);

	CHECK(geometry.term(0).factor(0, 2) == "11");
	CHECK(geometry.term(0).factor(0, 1) == "12");
	CHECK(geometry.term(1).factor(0, 2) == "21");
	CHECK(geometry.term(1).factor(0, 1) == "22");
}

TEST_CASE("unqualified GeometryFactor aliases geometry term zero direction zero",
          "[GeometryFactor]")
{
	Dmrg::InputCheck       inputCheck;
	InputNgType::Writeable writeable(inputCheck, oneTermInput("GeometryFactor=legacy;\n"));
	InputNgType::Readable  readable(writeable);
	GeometryType           geometry(readable);

	CHECK(geometry.term(0).factor(0, 2) == "legacy");
	CHECK(geometry.term(0).factor(0, 1).empty());
}

TEST_CASE("qualified GeometryFactor is accepted for a single geometry term", "[GeometryFactor]")
{
	Dmrg::InputCheck       inputCheck;
	InputNgType::Writeable writeable(inputCheck,
	                                 oneTermInput("gt0:dir0:GeometryFactor=explicit;\n"));
	InputNgType::Readable  readable(writeable);
	GeometryType           geometry(readable);

	CHECK(geometry.term(0).factor(0, 2) == "explicit");
	CHECK(geometry.term(0).factor(0, 1).empty());
}
