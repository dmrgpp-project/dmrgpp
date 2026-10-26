#include <PsimagLite/Ainur/Ainur.h>

#include <catch2/catch_test_macros.hpp>

TEST_CASE("Ainur expands macros in unused declarations", "[Ainur][macros]")
{
	PsimagLite::Ainur ainur("string! unusedMacro = \\macro;\n"
	                        "string macro = expanded;");
	std::string       value;
	ainur.readValue(value, "unusedMacro");

	CHECK(value == "expanded");
}
