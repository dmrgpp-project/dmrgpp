// SPDX-FileCopyrightText: Copyright (c) 2009-2016-2018, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [DMRG++ authors; see AUTHORS.]

/*! \file Su2Related.h
 *
 *
 */
#ifndef DMRG_SU2_RELATED_H
#define DMRG_SU2_RELATED_H
#include "InputCheck.h"
#include <PsimagLite/CrsMatrix.h>
#include <PsimagLite/InputNg.h>

namespace Dmrg {
// This is a structure, don't add member functions here!
struct Su2Related {
	Su2Related()
	    : offset(0) // setting to zero is necessary, because
	{ } // we always print offset
	// and when running Abelian
	// it might be undefined

	template <typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType& mres, SizeType, PsimagLite::String msg = "") const
	{
		return 0;
	}

	void read(PsimagLite::String label, PsimagLite::IoSerializer& ioSerializer)
	{
		ioSerializer.read(offset, label + "/offset");
		ioSerializer.read(source, label + "/source");
		ioSerializer.read(transpose, label + "/transpose");
	}

	void write(PsimagLite::String label, PsimagLite::IoSerializer& ioSerializer) const
	{
		ioSerializer.createGroup(label);
		ioSerializer.write(label + "/offset", offset);
		ioSerializer.write(label + "/source", source);
		ioSerializer.write(label + "/transpose", transpose);
	}

	SizeType                           offset;
	PsimagLite::Vector<SizeType>::Type source;
	PsimagLite::Vector<int>::Type      transpose;
};

std::istream& operator>>(std::istream&, Su2Related&);

std::ostream& operator<<(std::ostream&, const Su2Related&);

void send(Su2Related&, int, int, PsimagLite::MPI::CommType);

void recv(Su2Related&, int, int, PsimagLite::MPI::CommType);

void bcast(Su2Related&);

} // namespace Dmrg

#endif
