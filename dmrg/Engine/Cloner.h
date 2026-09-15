// SPDX-FileCopyrightText: Copyright (c) 2009, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 2.0.0]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file Cloner.h
 *
 *  A class to help repeating a task multiple times
 *
 */
#ifndef CLONER_HEADER_H
#define CLONER_HEADER_H

#include "String.h"
#include <PsimagLite/TypeToString.h>
#include <fstream>
#include <iostream>
#include <vector>

namespace Dmrg {
template <typename LineChangerType> class Cloner {
	static const SizeType LINE_LENGTH = 1024;

public:

	Cloner(const PsimagLite::String& infile,
	       const PsimagLite::String& outRoot,
	       const PsimagLite::String& ext)
	    : infile_(infile)
	    , outRoot_(outRoot)
	    , ext_(ext)
	{ }

	void push(const LineChangerType& lineChanger) { lineChanger_.push_back(lineChanger); }

	void createInputFile(SizeType ind) const
	{
		std::ifstream      fin(infile_.c_str());
		PsimagLite::String outfile = outRoot_ + ttos(ind) + ext_;
		std::ofstream      fout(outfile.c_str());
		char               line[LINE_LENGTH];
		while (!fin.eof()) {
			fin.getline(line, LINE_LENGTH);
			PsimagLite::String s(line);
			if (!procLine(s, ind))
				continue;
			fout << s << "\n";
		}
		fout.close();
		fin.close();
	}

private:

	bool procLine(PsimagLite::String& s, SizeType ind) const
	{
		for (SizeType i = 0; i < lineChanger_.size(); i++) {
			if (s.find(lineChanger_[i].string()) != PsimagLite::String::npos) {
				return lineChanger_[i].act(ind, s);
			}
		}
		return true;
	}

	PsimagLite::String                                 infile_, outRoot_, ext_;
	typename PsimagLite::Vector<LineChangerType>::Type lineChanger_;
}; // class Cloner

} // namespace Dmrg

/*@}*/
#endif // CLONER_HEADER_H
