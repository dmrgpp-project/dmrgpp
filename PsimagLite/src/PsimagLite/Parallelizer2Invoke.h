#ifndef PARALLELIZER2INVOKE_H
#define PARALLELIZER2INVOKE_H
#include <type_traits>

namespace PsimagLite {

template <typename LambdaType, typename ExecType>
constexpr bool parallelizer2WantsExec
    = std::is_invocable_v<const LambdaType&, SizeType, SizeType, const ExecType&>;

// Calls lambda(task, threadNum, exec) if the lambda accepts an execution
// space, else lambda(task, threadNum). Shared by all Parallelizer2 backends
// so that user lambdas have one signature regardless of backend.
template <typename LambdaType, typename ExecType>
void parallelizer2Invoke(const LambdaType& lambda,
                         SizeType          task,
                         SizeType          threadNum,
                         const ExecType&   exec)
{
	if constexpr (std::is_invocable_v<const LambdaType&, SizeType, SizeType, const ExecType&>)
		lambda(task, threadNum, exec);
	else
		lambda(task, threadNum);
}
} // namespace PsimagLite
#endif
