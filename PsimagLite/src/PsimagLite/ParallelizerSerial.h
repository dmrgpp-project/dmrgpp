#ifndef PARALLELIZER2SERIAL_H
#define PARALLELIZER2SERIAL_H
#include "CodeSectionParams.h"
#include "LoadBalancerDefault.h"
#include "ParallelizerInvoke.h"
#include "Vector.h"
#include <Kokkos_Core.hpp>

namespace PsimagLite {

template <typename = int> class Parallelizer {

public:

	Parallelizer(const CodeSectionParams& codeParams)
	{
		if (codeParams.npthreads != 1)
			throw RuntimeError("Please compile with -DUSE_PTHREADS\n");
	}

	template <typename SomeLambdaType>
	void parallelFor(SizeType start, SizeType end, const SomeLambdaType& lambda)
	{
		using ExecType = Kokkos::DefaultExecutionSpace;
		if constexpr (parallelizer2WantsExec<SomeLambdaType, ExecType>) {
			ExecType exec;
			for (SizeType i = start; i < end; ++i)
				parallelizer2Invoke(lambda, i, 0, exec);
		} else {
			for (SizeType i = start; i < end; ++i)
				lambda(i, 0);
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

	SizeType numberOfThreads() const { return 1; }

	String name() const { return "serial"; }
};
} // namespace PsimagLite
#endif // PARALLELIZER2SERIAL_H
