// SPDX-FileCopyrightText: Copyright (c) 2009-2016-2018, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [by G.A., Oak Ridge National Laboratory]

/*! \file Su2Related.cpp
 *
 *
 */
#include "Su2Related.h"

namespace Dmrg {

std::istream& operator>>(std::istream& is, Su2Related& x)
{
	is >> x.offset;
	return is;
}

std::ostream& operator<<(std::ostream& os, const Su2Related& x)
{
	os << x.offset << "\n";
	return os;
}

void send(Su2Related& su2Related, int root, int tag, PsimagLite::MPI::CommType mpiComm)
{
	PsimagLite::MPI::send(su2Related.offset, root, tag, mpiComm);
	PsimagLite::MPI::send(su2Related.source, root, tag + 1, mpiComm);
	PsimagLite::MPI::send(su2Related.transpose, root, tag + 2, mpiComm);
}

void recv(Su2Related& su2Related, int root, int tag, PsimagLite::MPI::CommType mpiComm)
{
	PsimagLite::MPI::recv(su2Related.offset, root, tag, mpiComm);
	PsimagLite::MPI::recv(su2Related.source, root, tag + 1, mpiComm);
	PsimagLite::MPI::recv(su2Related.transpose, root, tag + 2, mpiComm);
}

void bcast(Su2Related& su2Related)
{
	PsimagLite::MPI::bcast(su2Related.offset);
	PsimagLite::MPI::bcast(su2Related.source);
	PsimagLite::MPI::bcast(su2Related.transpose);
}

} // namespace Dmrg
