// SPDX-FileCopyrightText: Copyright (c) 2009-2020, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 2.]
// [see file AUTHORS for authorsL]

/** \ingroup DMRG */
/*@{*/

/*! \file InputCheckBase.h
 *
 *  InputChecking functions
 */
#ifndef PSI_INPUTCHECK_BASE_H
#define PSI_INPUTCHECK_BASE_H
#include <stdexcept>
#include <string>
#include <vector>

namespace PsimagLite {

class InputCheckBase {

	using VectorStringType = std::vector<std::string>;

public:

	std::string import() const { return ""; }

	bool check(const std::string&, const VectorStringType&, unsigned int) const
	{
		return false;
	}

	void check(const std::string&, const std::string&, unsigned int) { return; }

	bool checkSimpleLabel(const std::string&, unsigned int) const { return false; }

}; // class InputCheckBase
} // namespace PsimagLite

/*@}*/
#endif
