// SPDX-FileCopyrightText: Copyright (c) 2009-2013-2018, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file WaveFunctionTransfSu2.h
 *
 *  This class implements the wave function transformation, see PRL 77, 3633 (1996)
 *
 */

#ifndef WFT_BASE_H
#define WFT_BASE_H
#include "ProgramGlobals.h"
#include "WftOptions.h"
#include <PsimagLite/PackIndices.h>

namespace Dmrg {

template <typename DmrgWaveStructType_,
          typename VectorWithOffsetType_,
          typename OptionsType_,
          typename OneSiteSpacesType_>
class WaveFunctionTransfBase {

public:

	using DmrgWaveStructType     = DmrgWaveStructType_;
	using VectorWithOffsetType   = VectorWithOffsetType_;
	using PackIndicesType        = PsimagLite::PackIndices;
	using LeftRightSuperType     = typename DmrgWaveStructType::LeftRightSuperType;
	using BasisWithOperatorsType = typename DmrgWaveStructType::BasisWithOperatorsType;
	using BasisType              = typename BasisWithOperatorsType::BasisType;
	using VectorSizeType         = PsimagLite::Vector<SizeType>::Type;
	using WftOptionsType         = WftOptions<VectorWithOffsetType_, OptionsType_>;
	using OneSiteSpacesType      = OneSiteSpacesType_;

	virtual void transformVector(VectorWithOffsetType&       psiDest,
	                             const VectorWithOffsetType& psiSrc,
	                             const LeftRightSuperType&   lrs,
	                             const OneSiteSpacesType&    nk) const
	    = 0;

	virtual ~WaveFunctionTransfBase() { }
}; // class WaveFunctionTransfBase
} // namespace Dmrg

/*@}*/
#endif // WFT_BASE_H
