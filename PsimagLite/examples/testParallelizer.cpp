#include <PsimagLite/Concurrency.h>
#include <PsimagLite/LoadBalancerWeights.h>
#include <PsimagLite/Parallelizer.h>
#include <PsimagLite/Vector.h>

int main(int argc, char** argv)
{
	Kokkos::ScopeGuard scope_guard(argc, argv);

	constexpr unsigned int nthreads = 1;
	PsimagLite::Concurrency(&argc, &argv, nthreads);

	if (argc != 3) {
		std::cerr << "USAGE: " << argv[0] << " n n_threads\n";
		return 1;
	}

	const SizeType n       = atoi(argv[1]);
	const SizeType threads = atoi(argv[2]);

	PsimagLite::Vector<double>::Type v(n);

	PsimagLite::Parallelizer<> parallelizer(threads);

	std::cout << "Testing Parallelizer with " << parallelizer.name();
	std::cout << " and " << parallelizer.numberOfThreads() << " threads.\n";
	parallelizer.parallelFor(0,
	                         n,
	                         [&v](SizeType i, SizeType, const Kokkos::DefaultExecutionSpace&)
	                         { v[i] = i + 42; });

	/*
	for (SizeType i = 0; i < n; ++i) {
	        v[i] = i + 42;
	}
	*/

	std::cout << v[0] << " " << v[n - 1] << "\n";
}
