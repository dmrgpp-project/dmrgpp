#ifndef WFTHELPER_H
#define WFTHELPER_H
#include "OneSiteSpaces.hh"
#include "ParallelWftMany.h"
#include <OneSiteSpaces.hh>
#include <ParametersDmrgSolver.h>
#include <PsimagLite/Vector.h>

namespace Dmrg {

template <typename ModelType, typename VectorWithOffsetType, typename WaveFunctionTransfType>
class WftHelper {

public:

	using VectorVectorWithOffsetType = typename PsimagLite::Vector<VectorWithOffsetType*>::Type;
	using ModelHelperType            = typename ModelType::ModelHelperType;
	using LeftRightSuperType         = typename ModelHelperType::LeftRightSuperType;
	using VectorSizeType             = PsimagLite::Vector<SizeType>::Type;
	using OneSiteSpacesType          = OneSiteSpaces<ModelType>;
	using ParallelWftManyType        = ParallelWftMany<ModelType, WaveFunctionTransfType>;

	WftHelper(const ModelType&              model,
	          const LeftRightSuperType&     lrs,
	          const WaveFunctionTransfType& wft)
	    : model_(model)
	    , lrs_(lrs)
	    , wft_(wft)
	{ }

	void
	wftSome(VectorVectorWithOffsetType& tvs, SizeType site, SizeType begin, SizeType end) const
	{
		using ParallelizerType = PsimagLite::Parallelizer<ParallelWftManyType>;

		bool wants_parallel = wantsParallel();
		assert(end >= begin);
		SizeType total = end - begin;

		if (total == 0)
			return; // <-- EARLY EXIT HERE

		SizeType threads = (wants_parallel)
		    ? std::min(total, PsimagLite::Concurrency::codeSectionParams.npthreads)
		    : 1;

		PsimagLite::CodeSectionParams codeSectionParams(threads);
		ParallelizerType              threadedCtor(codeSectionParams);
		std::cout << "HERE " << threads << " threads\n";
		std::vector<SizeType> index_map(total);
		for (SizeType i = 0; i < total; ++i) {
			index_map[i] = i + begin;
		}

		ProgramGlobals::DirectionEnum dir
		    = ProgramGlobals::DirectionEnum::EXPAND_SYSTEM; // FIXME!
		OneSiteSpacesType   oneSiteSpaces(site, dir, model_);
		ParallelWftManyType helper(tvs, oneSiteSpaces, index_map, wft_, lrs_);

		threadedCtor.loopCreate(helper);
	}

	void wftOneVector(VectorWithOffsetType&       phiNew,
	                  const VectorWithOffsetType& src,
	                  SizeType                    site) const
	{
		ProgramGlobals::DirectionEnum dir
		    = ProgramGlobals::DirectionEnum::EXPAND_SYSTEM; // FIXME!
		OneSiteSpacesType one_site_space(site, dir, model_);
		ParallelWftManyType::wftOneVector(phiNew, src, one_site_space, wft_, lrs_);
	}

private:

	bool wantsParallel() const
	{
		const auto& opts = model_.params().options;
		if (opts.isSet("parallelwftmany")) {
			assert(opts.isSet("nowft") || opts.isSet("wftNoAccel"));
			return true;
		}

		return false;
	}

	const ModelType&              model_;
	const LeftRightSuperType&     lrs_;
	const WaveFunctionTransfType& wft_;
};
}
#endif // WFTHELPER_H
