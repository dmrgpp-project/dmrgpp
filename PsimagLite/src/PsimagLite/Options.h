// SPDX-FileCopyrightText: Copyright (c) 2012 , UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup PsimagLite */
/*@{*/

/*! \file Options.h
 *
 *  Options lines
 */

#ifndef OPTIONS_HEADER_H
#define OPTIONS_HEADER_H

#include "PsimagLite.h"
#include "Vector.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <stdexcept>

namespace PsimagLite {

class Options {

public:

	class Writeable {

		using VectorStringType = Vector<String>::Type;

	public:

		enum
		{
			DISABLED,
			PERMISSIVE,
			STRICT
		};

		Writeable(VectorStringType& registeredOptions, SizeType mode)
		    : registeredOptions_(registeredOptions)
		    , mode_(mode)
		{ }

		void set(Vector<String>::Type& optsThatAreSet, String opts)
		{
			if (mode_ == DISABLED)
				return;
			split(optsThatAreSet, opts, ",");
			for (SizeType i = 0; i < optsThatAreSet.size(); i++) {
				bool b = (find(registeredOptions_.begin(),
				               registeredOptions_.end(),
				               optsThatAreSet[i])
				          == registeredOptions_.end());
				if (!b)
					continue;

				String s(__FILE__);
				s += ": Unknown option " + optsThatAreSet[i] + "\n";
				if (mode_ == PERMISSIVE)
					std::cerr << " *** WARNING **: " << s;
				if (mode_ == STRICT)
					throw RuntimeError(s.c_str());
			}
		}

	private:

		Vector<String>::Type registeredOptions_;
		SizeType             mode_;
	}; // class Writeable

	class Readable {

	public:

		Readable(Writeable& optsWrite, const String& optsString)
		{
			optsWrite.set(optsThatAreSet_, optsString);
		}

		bool isSet(const String& thisOption) const
		{
			bool b = (find(optsThatAreSet_.begin(), optsThatAreSet_.end(), thisOption)
			          == optsThatAreSet_.end());
			return (!b);
		}

	private:

		Vector<String>::Type optsThatAreSet_;

	}; // class Readable

}; // class OptionsWriteable
} // namespace PsimagLite
/*@}*/
#endif // OPTIONS_HEADER_H
