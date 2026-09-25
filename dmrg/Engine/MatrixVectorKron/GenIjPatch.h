// SPDX-FileCopyrightText: Copyright (c) 2012, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 2.0.0]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file GenIjPatch.h
 *
 *
 */

#ifndef GEN_IJ_PATCH_HEADER_H
#define GEN_IJ_PATCH_HEADER_H

#include <PsimagLite/Vector.h>
#include <cassert>

namespace Dmrg {

template <typename LeftRightSuperType_> class GenIjPatch {

public:

	using LeftRightSuperType = LeftRightSuperType_;
	using BasisType          = typename LeftRightSuperType::BasisType;
	using QnType             = typename BasisType::QnType;
	using VectorSizeType     = PsimagLite::Vector<SizeType>::Type;

	enum LeftOrRightEnumType
	{
		LEFT  = 0,
		RIGHT = 1
	};

	GenIjPatch(const LeftRightSuperType& lrs, const QnType& target)
	    : lrs_(lrs)
	    , qn_(target)
	{
		for (SizeType i = 0; i < lrs.left().partition() - 1; i++) {
			for (SizeType j = 0; j < lrs.right().partition() - 1; ++j) {

				if (QnType(lrs.left().qnEx(i), lrs.right().qnEx(j)) != target)
					continue;

				patchesLeft_.push_back(i);
				patchesRight_.push_back(j);
			}
		}
	}

	const QnType& qn() const { return qn_; }

	const VectorSizeType& operator()(LeftOrRightEnumType leftOrRight) const
	{
		return (leftOrRight == LEFT) ? patchesLeft_ : patchesRight_;
	}

	const LeftRightSuperType& lrs() const { return lrs_; }

private:

	const LeftRightSuperType& lrs_;
	const QnType&             qn_;
	VectorSizeType            patchesLeft_;
	VectorSizeType            patchesRight_;

}; // class GenIjPatch
} // namespace PsimagLite

/*@}*/

#endif // GEN_IJ_PATCH_HEADER_H
