// SPDX-FileCopyrightText: Copyright (c) 2014, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup PsimagLite */
/*@{*/

/*! \file Star.h
 *
 *  DOC NEEDED FIXME
 */
#ifndef GEOMETRY_STAR_H
#define GEOMETRY_STAR_H
#include "GeometryBase.h"

namespace PsimagLite {

template <typename ComplexOrRealType, typename InputType>
class Star : public GeometryBase<ComplexOrRealType, InputType> {

public:

	enum
	{
		DIRECTION_S
	};

	static const SizeType CENTER = 0;

	Star() { }

	Star(SizeType linSize, InputType&)
	    : linSize_(linSize)
	{ }

	SizeType maxConnections() const override { return linSize_ - 1; }

	SizeType dirs() const override { return 1; }

	SizeType handle(SizeType i, SizeType j) const override
	{
		assert(i == CENTER || j == CENTER);
		SizeType k      = (i == CENTER) ? j : i;
		SizeType offset = (k < CENTER) ? 0 : 1;
		assert(k >= offset);
		return k - offset;
	}

	SizeType getVectorSize(SizeType dirId) const override
	{
		if (dirId != DIRECTION_S)
			throw RuntimeError("Star must have direction 0\n");
		return linSize_ - 1;
	}

	bool connected(SizeType i1, SizeType i2) const override
	{
		if (i1 == i2)
			return false;
		return (i1 == CENTER || i2 == CENTER);
	}

	// assumes i1 and i2 are connected
	SizeType calcDir(SizeType, SizeType) const override { return DIRECTION_S; }

	bool fringe(SizeType i, SizeType, SizeType emin) const override { return (i == emin); }

	// siteNew2 is fringe in the environment
	SizeType getSubstituteSite(SizeType, SizeType, SizeType siteNew2) const override
	{
		return siteNew2;
	}

	String label() const override { return "star"; }

	SizeType findReflection(SizeType) const override { throw RuntimeError("findReflection\n"); }

	SizeType length(SizeType i) const override
	{
		assert(i == 0);
		return linSize_;
	}

	SizeType translate(SizeType, SizeType dir, SizeType) const override
	{
		assert(dir == 0);

		throw RuntimeError("translate\n");
	}

private:

	SizeType linSize_;
}; // class Star
} // namespace PsimagLite

/*@}*/
#endif // GEOMETRY_STAR_H
