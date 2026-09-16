// SPDX-FileCopyrightText: Copyright (c) 2009, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 2.0.0]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file LineChangerLinear.h
 *
 *
 *
 */
#ifndef LINE_CHANGER_LINEAR_H
#define LINE_CHANGER_LINEAR_H

#include "String.h"
#include <PsimagLite/TypeToString.h>
#include <vector>

namespace Dmrg {
template <typename ValueType> class LineChangerLinear {
public:

	LineChangerLinear(const PsimagLite::String& match,
	                  ValueType                 step,
	                  ValueType                 init,
	                  const PsimagLite::String& pre,
	                  const PsimagLite::String& post)
	    : match_(match)
	    , step_(step)
	    , init_(init)
	    , pre_(pre)
	    , post_(post)
	{ }

	const PsimagLite::String& string() const { return match_; }

	bool act(SizeType i, PsimagLite::String& line) const
	{
		ValueType val = i * step_ + init_;
		line          = match_ + pre_ + ttos(val) + post_;
		return true;
	}

private:

	PsimagLite::String match_;
	ValueType          step_, init_;
	PsimagLite::String pre_, post_;
}; // class LineChangerLinear

} // namespace Dmrg

/*@}*/
#endif // LINE_CHANGER_LINEAR_H
