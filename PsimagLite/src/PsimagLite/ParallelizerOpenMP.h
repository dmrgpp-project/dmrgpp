#ifndef PARALLELIZER2OPENMP_H
#define PARALLELIZER2OPENMP_H
#include "CodeSectionParams.h"
#include "LoadBalancerDefault.h"
#include "ParallelizerInvoke.h"
#include "Vector.h"
#include <Kokkos_Core.hpp>
#include <omp.h>

namespace PsimagLite {

template <typename = int> class Parallelizer {

public:

	Parallelizer(const CodeSectionParams& codeParams)
	    : threads_(codeParams.npthreads)
	{
		omp_set_num_threads(codeParams.npthreads);
	}

	template <typename SomeLambdaType>
	void parallelFor(SizeType start, SizeType end, const SomeLambdaType& lambda)
	{
		using ExecType = Kokkos::DefaultExecutionSpace;
		if constexpr (parallelizer2WantsExec<SomeLambdaType, ExecType>) {
			ExecType exec;
#pragma omp parallel for
			for (SizeType i = start; i < end; ++i)
				parallelizer2Invoke(lambda, i, omp_get_thread_num(), exec);
		} else {
#pragma omp parallel for
			for (SizeType i = start; i < end; ++i)
				lambda(i, omp_get_thread_num());
		}
	}

	// weights and balancers only matter for the pthread backend; accepted
	// here so that callers compile unchanged with any backend
	template <typename SomeLambdaType>
	void parallelFor(SizeType              start,
	                 SizeType              end,
	                 const SomeLambdaType& lambda,
	                 const LoadBalancerDefault::VectorSizeType&)
	{
		parallelFor(start, end, lambda);
	}

	template <typename SomeLambdaType, typename SomeLoadBalancer>
	void parallelFor(SizeType              start,
	                 SizeType              end,
	                 const SomeLambdaType& lambda,
	                 const SomeLoadBalancer&)
	{
		parallelFor(start, end, lambda);
	}

	SizeType numberOfThreads() const { return omp_get_num_threads(); }

	String name() const { return "openmp"; }

private:

	SizeType threads_;
};
} // namespace PsimagLite
#endif // PARALLELIZER2OPENMP_H
