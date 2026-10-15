/*
Copyright (c) 2009-2015, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 5.]
[by G.A., Oak Ridge National Laboratory]

UT Battelle Open Source Software License 11242008

OPEN SOURCE LICENSE

Subject to the conditions of this License, each
contributor to this software hereby grants, free of
charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), a
perpetual, worldwide, non-exclusive, no-charge,
royalty-free, irrevocable copyright license to use, copy,
modify, merge, publish, distribute, and/or sublicense
copies of the Software.

1. Redistributions of Software must retain the above
copyright and license notices, this list of conditions,
and the following disclaimer.  Changes or modifications
to, or derivative works of, the Software should be noted
with comments and the contributor and organization's
name.

2. Neither the names of UT-Battelle, LLC or the
Department of Energy nor the names of the Software
contributors may be used to endorse or promote products
derived from this software without specific prior written
permission of UT-Battelle.

3. The software and the end-user documentation included
with the redistribution, with or without modification,
must include the following acknowledgment:

"This product includes software produced by UT-Battelle,
LLC under Contract No. DE-AC05-00OR22725  with the
Department of Energy."

*********************************************************
DISCLAIMER

THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT OWNER, CONTRIBUTORS, UNITED STATES GOVERNMENT,
OR THE UNITED STATES DEPARTMENT OF ENERGY BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED
STATES DEPARTMENT OF ENERGY, NOR THE COPYRIGHT OWNER, NOR
ANY OF THEIR EMPLOYEES, REPRESENTS THAT THE USE OF ANY
INFORMATION, DATA, APPARATUS, PRODUCT, OR PROCESS
DISCLOSED WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*********************************************************

 */
/** \ingroup DMRG */
/*@{*/
/** \file ParallelWftMany.h
 */

#ifndef DMRG_PARALLEL_WFT_MANY_H
#define DMRG_PARALLEL_WFT_MANY_H

#include <OneSiteSpaces.hh>
#include <PsimagLite/Concurrency.h>
#include <vector>

namespace Dmrg {

template <typename ModelType, typename WaveFunctionTransfType> class ParallelWftMany {

	using ComplexOrRealType          = typename ModelType::ComplexOrRealType;
	using ConcurrencyType            = PsimagLite::Concurrency;
	using VectorWithOffsetType       = VectorWithOffsets<ComplexOrRealType>;
	using VectorVectorWithOffsetType = std::vector<VectorWithOffsetType*>;
	using OneSiteSpacesType          = OneSiteSpaces<ModelType>;
	using ModelHelperType            = typename ModelType::ModelHelperType;
	using LeftRightSuperType         = typename ModelHelperType::LeftRightSuperType;

public:

	ParallelWftMany(VectorVectorWithOffsetType&   targetVectors,
	                const OneSiteSpacesType&      one_site_spaces,
	                const std::vector<SizeType>&  index_map,
	                const WaveFunctionTransfType& wft,
	                const LeftRightSuperType&     lrs)
	    : targetVectors_(targetVectors)
	    , one_site_spaces_(one_site_spaces)
	    , index_map_(index_map)
	    , wft_(wft)
	    , lrs_(lrs)
	{ }

	void doTask(SizeType ix, SizeType /* threadNum */)

	{
		assert(ix < index_map_.size());
		SizeType index = index_map_[ix];

		assert(index < targetVectors_.size());
		const VectorWithOffsetType& src = *targetVectors_[index];
		if (src.size() == 0)
			return; // <-- EARLY EXIT HERE

		VectorWithOffsetType phiNew;
		wftOneVector(phiNew, src, one_site_spaces_, wft_, lrs_);
		*targetVectors_[index] = phiNew;
		std::cout << "HERE: task " << ix << " out of " << tasks() << "\n";
	}

	SizeType tasks() const { return index_map_.size(); }

	static void wftOneVector(VectorWithOffsetType&         phiNew,
	                         const VectorWithOffsetType&   src,
	                         const OneSiteSpacesType&      one_site_spaces,
	                         const WaveFunctionTransfType& wft,
	                         const LeftRightSuperType&     lrs)
	{
		phiNew.populateFromQns(src, lrs.super());

		// OK, now that we got the partition number right, let's wft:
		wft.setInitialVector(phiNew, src, lrs, one_site_spaces);
	}

private:

	VectorVectorWithOffsetType&   targetVectors_;
	const OneSiteSpacesType&      one_site_spaces_;
	std::vector<SizeType>         index_map_;
	const WaveFunctionTransfType& wft_;
	const LeftRightSuperType&     lrs_;
}; // class ParallelWftMany
} // namespace Dmrg

/*@}*/
#endif // DMRG_PARALLEL_WFT_MANY_H
