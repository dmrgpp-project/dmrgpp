// To be included by Mpi.cpp ONLY
#include "MpiNo.h"

namespace PsimagLite {

namespace MPI {

	int comm_world = 0;
	int sum        = 0;

	void init(int*, char***) { }

	void finalize() { }

	SizeType commSize(CommType) { return 1; }

	SizeType commRank(CommType) { return 0; }

	int barrier(CommType) { return 0; }

	bool hasMpi() { return false; }

	void info(std::ostream&) { }

	void version(std::ostream&) { }

} // namespace MPI

} // namespace PsimagLite
