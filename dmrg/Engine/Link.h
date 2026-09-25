// SPDX-FileCopyrightText: Copyright (c) 2009-2013, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 2.0.0]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file Link.h
 *
 *  DOC NEEDED FIXME
 */
#ifndef LINK_H
#define LINK_H

#include "MetaOpForConnection.hh"
#include "ProgramGlobals.h"

namespace Dmrg {

template <typename FieldType> struct Link {

	using PairSizeType               = std::pair<SizeType, SizeType>;
	using RealType                   = typename PsimagLite::Real<FieldType>::Type;
	using PairMetaOForConnectoinType = std::pair<MetaOpForConnection, MetaOpForConnection>;

	Link(PairMetaOForConnectoinType         finalIndices1,
	     ProgramGlobals::ConnectionEnum     type1,
	     const FieldType&                   value1,
	     ProgramGlobals::FermionOrBosonEnum fOb,
	     SizeType                           aM,
	     RealType                           aF,
	     SizeType                           cat)
	    : pairMetaOps(finalIndices1)
	    , type(type1)
	    , value(value1)
	    , fermionOrBoson(fOb)
	    , angularMomentum(aM)
	    , angularFactor(aF)
	    , category(cat)
	{ }

	PairMetaOForConnectoinType         pairMetaOps;
	ProgramGlobals::ConnectionEnum     type;
	FieldType                          value;
	ProgramGlobals::FermionOrBosonEnum fermionOrBoson;
	SizeType                           angularMomentum;
	RealType                           angularFactor;
	SizeType                           category;
}; // struct Link
} // namespace Dmrg

/*@}*/
#endif // LINK_H
