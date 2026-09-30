#include "ProcOmegas.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <cstdio>
#include <fstream>

namespace {

struct OmegaParamsStub { };

using Parser     = Dmrg::ProcOmegas<double, OmegaParamsStub>;
using VectorType = Parser::VectorRealType;

void writeLog(const char* filename, const char* contents)
{
	std::ofstream fout(filename);
	REQUIRE(fout.good());
	fout << contents;
}

} // namespace

TEST_CASE("ProcOmegas reads complete P2 and P3 data", "[ProcOmegas]")
{
	const char* filename = "test_ProcOmegas_complete.cout";
	writeLog(filename,
	         "0 1.5 0.0 <gs|c'|P2> 0.0\n"
	         "1 2.5 0.0 <gs|c'|P2> 0.0\n"
	         "0 3.5 0.0 <gs|c'|P3> 0.0\n"
	         "1 4.5 0.0 <gs|c'|P3> 0.0\n");

	VectorType p2({ -1.0, -1.0 });
	VectorType p3({ -1.0, -1.0 });
	Parser::correctionVectorRead(p2, p3, filename);

	CHECK(p2[0] == 1.5);
	CHECK(p2[1] == 2.5);
	CHECK(p3[0] == 3.5);
	CHECK(p3[1] == 4.5);
	std::remove(filename);
}

TEST_CASE("ProcOmegas rejects incomplete P2 data", "[ProcOmegas]")
{
	const char* filename = "test_ProcOmegas_missing-p2.cout";
	writeLog(filename,
	         "0 1.5 0.0 <gs|c'|P2> 0.0\n"
	         "0 3.5 0.0 <gs|c'|P3> 0.0\n"
	         "1 4.5 0.0 <gs|c'|P3> 0.0\n");

	VectorType p2(2);
	VectorType p3(2);
	CHECK_THROWS_WITH(Parser::correctionVectorRead(p2, p3, filename),
	                  "Undefined P2 value for site= 1 file= test_ProcOmegas_missing-p2.cout\n");
	std::remove(filename);
}

TEST_CASE("ProcOmegas rejects incomplete P3 data", "[ProcOmegas]")
{
	const char* filename = "test_ProcOmegas_missing-p3.cout";
	writeLog(filename,
	         "0 1.5 0.0 <gs|c'|P2> 0.0\n"
	         "1 2.5 0.0 <gs|c'|P2> 0.0\n"
	         "0 3.5 0.0 <gs|c'|P3> 0.0\n");

	VectorType p2(2);
	VectorType p3(2);
	CHECK_THROWS_WITH(Parser::correctionVectorRead(p2, p3, filename),
	                  "Undefined P3 value for site= 1 file= test_ProcOmegas_missing-p3.cout\n");
	std::remove(filename);
}
