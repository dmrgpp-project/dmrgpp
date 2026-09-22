#include <PsimagLite/TextAndNumbersChecker.hpp>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

enum class ComparisonMode
{
	GENERIC,
	IGNORE_LINE_PREFIX,
	OBSERVABLE
};

struct ParsedOption {
	ComparisonMode mode;
	std::string    value;
};

ParsedOption parseOption(const char* text)
{
	const std::string argument(text);
	const std::string ignored_line_prefix_option = "--ignore-line-prefix=";
	const std::string insitu_label_option        = "--insitu-label=";

	if (argument.compare(0, ignored_line_prefix_option.size(), ignored_line_prefix_option)
	    == 0) {
		if (argument.size() == ignored_line_prefix_option.size())
			throw std::invalid_argument("Invalid option: \"" + argument + "\"");
		return { ComparisonMode::IGNORE_LINE_PREFIX,
			 argument.substr(ignored_line_prefix_option.size()) };
	}

	if (argument.compare(0, insitu_label_option.size(), insitu_label_option) == 0) {
		if (argument.size() == insitu_label_option.size())
			throw std::invalid_argument("Invalid option: \"" + argument + "\"");
		return { ComparisonMode::OBSERVABLE, argument.substr(insitu_label_option.size()) };
	}

	throw std::invalid_argument("Invalid option: \"" + argument + "\"");
}

double parseTolerance(const char* text)
{
	const std::string token(text);
	std::size_t       parsed = 0;
	try {
		const double value = std::stod(token, &parsed);
		if (parsed != token.size() || !std::isfinite(value))
			throw std::invalid_argument("not a finite real number");
		return value;
	} catch (const std::invalid_argument&) {
		throw std::invalid_argument("Invalid tolerance: \"" + token + "\"");
	} catch (const std::out_of_range&) {
		throw std::invalid_argument("Tolerance is out of range: \"" + token + "\"");
	}
}

} // namespace

int main(int argc, char* argv[])
{
	if (argc != 4 && argc != 5) {
		std::cerr << "Usage: " << argv[0]
		          << " FILE1 FILE2 TOLERANCE "
		             "[--ignore-line-prefix=PREFIX|--insitu-label=LABEL]\n";
		return EXIT_FAILURE;
	}

	try {
		const ParsedOption option = (argc == 5)
		    ? parseOption(argv[4])
		    : ParsedOption { ComparisonMode::GENERIC, "" };
		const std::string  ignored_line_prefix
		    = (option.mode == ComparisonMode::IGNORE_LINE_PREFIX) ? option.value : "";
		const PsimagLite::TextAndNumbersChecker checker(parseTolerance(argv[3]),
		                                                ignored_line_prefix);

		if (option.mode == ComparisonMode::OBSERVABLE)
			checker.runObservables(argv[1], argv[2], option.value);
		else
			checker.run(argv[1], argv[2]);
	} catch (const std::exception& error) {
		std::cerr << error.what() << '\n';
		return EXIT_FAILURE;
	}
}
