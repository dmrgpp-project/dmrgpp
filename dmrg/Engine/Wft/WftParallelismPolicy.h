#ifndef WFT_PARALLELISM_POLICY_H
#define WFT_PARALLELISM_POLICY_H
#include <PsimagLite/Concurrency.h>
#include <PsimagLite/TypeToString.h>

namespace Dmrg {

class WftParallelismPolicy {

public:

	template <typename WftOptionsType>
	static void check(typename WftOptionsType::AccelEnum   accel,
	                  const PsimagLite::CodeSectionParams& outer,
	                  const PsimagLite::CodeSectionParams& blocks,
	                  const PsimagLite::CodeSectionParams& patchesCtor,
	                  const PsimagLite::CodeSectionParams& patchesTransform)
	{
		if (outer.npthreads == 1)
			return;

		if (accel == WftOptionsType::ACCEL_BLOCKS && blocks.npthreads > 1)
			throw PsimagLite::RuntimeError(message("blocks", outer, blocks));

		if (accel != WftOptionsType::ACCEL_PATCHES)
			return;

		const bool constructorConflict = (patchesCtor.npthreads > 1);
		const bool transformConflict   = (patchesTransform.npthreads > 1);
		if (!constructorConflict && !transformConflict)
			return;

		PsimagLite::String msg = prefix("patches", outer);
		if (constructorConflict)
			msg += patchesCtor.name + "=" + ttos(patchesCtor.npthreads);
		if (constructorConflict && transformConflict)
			msg += ", ";
		if (transformConflict)
			msg += patchesTransform.name + "=" + ttos(patchesTransform.npthreads);
		msg += "; set one of these thread counts to 1\n";
		throw PsimagLite::RuntimeError(msg);
	}

private:

	static PsimagLite::String prefix(const PsimagLite::String&            mode,
	                                 const PsimagLite::CodeSectionParams& outer)
	{
		return outer.name + "=" + ttos(outer.npthreads) + " conflicts with accelerated WFT "
		    + mode + " section ";
	}

	static PsimagLite::String message(const PsimagLite::String&            mode,
	                                  const PsimagLite::CodeSectionParams& outer,
	                                  const PsimagLite::CodeSectionParams& inner)
	{
		return prefix(mode, outer) + inner.name + "=" + ttos(inner.npthreads)
		    + "; set one of these thread counts to 1\n";
	}
};
} // namespace Dmrg
#endif // WFT_PARALLELISM_POLICY_H
