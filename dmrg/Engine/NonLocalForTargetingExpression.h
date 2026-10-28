#ifndef NONLOCALFORTARGETINGEXPRESSION_H
#define NONLOCALFORTARGETINGEXPRESSION_H
#include "AuxForTargetingExpression.h"
#include "TargetParamsBase.h"
#include <PsimagLite/OneOperatorSpec.h>
#include <PsimagLite/ProgressIndicator.h>
#include <string>
#include <vector>

namespace Dmrg {

template <typename TargetingBaseType> class NonLocalForTargetingExpression {

public:

	using AuxiliaryType               = AuxForTargetingExpression<TargetingBaseType>;
	using VectorWithOffsetType        = typename TargetingBaseType::VectorWithOffsetType;
	using ModelType                   = typename TargetingBaseType::ModelType;
	using TargetParamsType            = typename TargetingBaseType::TargetParamsType;
	using RealType                    = typename TargetParamsType::RealType;
	using ApplyOperatorExpressionType = typename TargetingBaseType::ApplyOperatorExpressionType;
	using GroupOfOneTimeEvolutionsType = typename AuxiliaryType::GroupOfOneTimeEvolutionsType;
	using OneTimeEvolutionType = typename GroupOfOneTimeEvolutionsType::OneTimeEvolutionType;
	using OneOperatorSpecType  = PsimagLite::OneOperatorSpec;
	using SiteSplitType        = typename OneOperatorSpecType::SiteSplit;
	using TimeVectorsBaseType  = typename ApplyOperatorExpressionType::TimeVectorsBaseType;

	NonLocalForTargetingExpression(const AuxiliaryType& aux)
	    : aux_(aux)
	    , progress_("NonLocalForTargetingExpression")
	{ }

	bool timeEvolve(const std::string&   name,
	                const SiteSplitType& siteSplit,
	                const std::string&   srcKet,
	                SizeType             site)
	{
		if (siteSplit.hasSiteString)
			err("Global operators cannot have a site\n");

		const SizeType firstIndex = aux_.pIndexOutput();
		if (firstIndex >= aux_.pVectors().origPvectors())
			err("Cannot set an NGST to a temporary vector\n");

		OneTimeEvolutionType* oneTimeEvolution
		    = aux_.timeEvolve().findThisEvolution(firstIndex);
		const bool                  isNewEvolution = (oneTimeEvolution == nullptr);
		const PsimagLite::String    legacySrcKet(srcKet.c_str());
		const VectorWithOffsetType* srcVwo = &aux_.getCurrentVectorConst(legacySrcKet);
		if (srcVwo->size() == 0 && !oneTimeEvolution)
			return false;

		AuxiliaryType* auxPtr = const_cast<AuxiliaryType*>(&aux_);
		if (isNewEvolution) {
			const SizeType nsites
			    = aux_.pVectors().aoe().model().superGeometry().numberOfSites();
			oneTimeEvolution = new OneTimeEvolutionType(
			    srcKet, name, nsites, aux_.Eg(), aux_.getCurrentTime(legacySrcKet));
		}

		using LegacyVectorSizeType = typename ModelType::VectorSizeType;
		using LegacyVectorRealType = typename PsimagLite::Vector<RealType>::Type;
		LegacyVectorRealType chebyTransform(oneTimeEvolution->chebyTransform().begin(),
		                                    oneTimeEvolution->chebyTransform().end());
		auxPtr->pVectors().initTimeVectors(
		    oneTimeEvolution->timeSteps(),
		    oneTimeEvolution->tau(),
		    PsimagLite::String(oneTimeEvolution->algorithm().c_str()),
		    chebyTransform);
		if (isNewEvolution) {
			oneTimeEvolution->initialize(firstIndex, *srcVwo, aux_.pVectors());
			aux_.timeEvolve().pushBack(oneTimeEvolution);
		}

		const bool timeHasAdvanced = advanceInTimeOrNot(*oneTimeEvolution,
		                                                oneTimeEvolution->advanceEach(),
		                                                oneTimeEvolution->depends(),
		                                                site,
		                                                oneTimeEvolution->tau());

		assert(oneTimeEvolution->indices().size() > 1);
		const SizeType              last         = oneTimeEvolution->indices().size() - 1;
		const SizeType              advanceOrNot = (timeHasAdvanced) ? last : 0;
		const SizeType              firstOrLast = oneTimeEvolution->indices()[advanceOrNot];
		const VectorWithOffsetType* phi         = isNewEvolution
		            ? srcVwo
		            : new VectorWithOffsetType(aux_.pVectors().aoe().targetVectors(firstOrLast));

		const LegacyVectorSizeType indices(oneTimeEvolution->indices().begin(),
		                                   oneTimeEvolution->indices().end());
		const LegacyVectorSizeType block1(1, site);
		static const bool          allOperatorsApplied = true;
		static const bool          isLastCall          = true;
		auxPtr->pVectors().aoeNonConst().calcTimeVectors(indices,
		                                                 oneTimeEvolution->energyForExp(),
		                                                 *phi,
		                                                 aux_.direction(),
		                                                 allOperatorsApplied,
		                                                 false,
		                                                 block1,
		                                                 isLastCall);

		assert(oneTimeEvolution->indices().size() == oneTimeEvolution->timeSteps());
		const RealType baseTime = oneTimeEvolution->time();
		for (SizeType i = 0; i < oneTimeEvolution->indices().size(); ++i) {
			assert(oneTimeEvolution->timeSteps() > 1);
			const RealType offset
			    = i * oneTimeEvolution->tau() / (oneTimeEvolution->timeSteps() - 1);
			auxPtr->pVectors().setTime(oneTimeEvolution->indices()[i],
			                           baseTime + offset);
		}

		if (!isNewEvolution)
			delete phi;

		return true;
	}

private:

	bool advanceInTimeOrNot(OneTimeEvolutionType&        oneTimeEvolution,
	                        SizeType                     advanceEach,
	                        const std::vector<SizeType>& depends,
	                        SizeType                     site,
	                        RealType                     tau)
	{
		if (!passDepends(depends))
			return false;

		static const bool advanceOnlyAtBorder = true;
		const SizeType    sites
		    = aux_.pVectors().aoe().model().superGeometry().numberOfSites();
		const bool weAreAtBorder   = (site == 0 || site == sites - 1);
		const bool dontAdvance     = (advanceOnlyAtBorder && !weAreAtBorder);
		bool       timeHasAdvanced = false;
		if (advanceEach > 0 && oneTimeEvolution.timesWithoutAdvancement() >= advanceEach
		    && !dontAdvance) {
			oneTimeEvolution.resetTimesWithoutAdvancement();
			oneTimeEvolution.advanceTime(tau);
			timeHasAdvanced           = true;
			const int sourceToDestroy = oneTimeEvolution.sourceToDestroy();
			if (sourceToDestroy >= 0)
				aux_.pVectors().aoeNonConst().destroyPvector(sourceToDestroy);
		} else {
			oneTimeEvolution.incrementTimesWithoutAdvancement();
		}

		TimeVectorsBaseType* ptr
		    = const_cast<TimeVectorsBaseType*>(&aux_.pVectors().aoe().timeVectors());
		ptr->setCurrentTime(oneTimeEvolution.time());

		PsimagLite::OstringStream msgg2(std::cout.precision());
		msgg2() << "Steps without advance: " << oneTimeEvolution.timesWithoutAdvancement();
		msgg2() << " site=" << site << " currenTime=" << oneTimeEvolution.time();
		if (oneTimeEvolution.timesWithoutAdvancement() > 0)
			progress_.printline(msgg2, std::cout);
		return timeHasAdvanced;
	}

	bool passDepends(const std::vector<SizeType>& depends)
	{
		for (SizeType index : depends)
			if (!aux_.pVectors()(index).isDone())
				return false;

		return true;
	}

	const AuxiliaryType&          aux_;
	PsimagLite::ProgressIndicator progress_;
};
}
#endif // NONLOCALFORTARGETINGEXPRESSION_H