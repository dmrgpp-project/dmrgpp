#include <PsimagLite/Io/IoNg.h>
#include <PsimagLite/PsimagLite.h>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char* argv[])
{
	if (argc < 5) {
		std::cerr << "Usage: " << argv[0]
		          << " file dataset tolerance expected [expected ...]\n";
		return 1;
	}

	try {
		const double tolerance = std::stod(argv[3]);
		if (tolerance < 0) {
			std::cerr << "Tolerance must be nonnegative\n";
			return 1;
		}

		std::vector<double> expected;
		expected.reserve(argc - 4);
		for (int i = 4; i < argc; ++i)
			expected.push_back(std::stod(argv[i]));

		PsimagLite::IoNg::In io(argv[1]);
		std::vector<double>  actual;
		io.read(actual, argv[2]);

		if (actual.size() != expected.size()) {
			std::cerr << "Dataset " << argv[2] << " has " << actual.size()
			          << " values; expected " << expected.size() << "\n";
			return 1;
		}

		for (SizeType i = 0; i < actual.size(); ++i) {
			if (std::abs(actual[i] - expected[i]) <= tolerance)
				continue;

			std::cerr << "Dataset " << argv[2] << " value " << i << " is " << actual[i]
			          << "; expected " << expected[i] << " within " << tolerance
			          << "\n";
			return 1;
		}
	} catch (const std::exception& error) {
		std::cerr << "check_hdf5_vector: " << error.what() << "\n";
		return 1;
	} catch (...) {
		std::cerr << "check_hdf5_vector: HDF5 read failed\n";
		return 1;
	}

	return 0;
}
