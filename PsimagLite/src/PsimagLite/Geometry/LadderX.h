// SPDX-FileCopyrightText: Copyright (c) 2009-2013-2018, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 2.]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup PsimagLite */
/*@{*/

/*! \file LadderX.h
 *
 *  DOC NEEDED FIXME
 */
#ifndef LADDERX_H
#define LADDERX_H

#include "Ladder.h"
#include <stdexcept>

namespace PsimagLite {

template <typename ComplexOrRealType, typename InputType>
class LadderX : public GeometryBase<ComplexOrRealType, InputType> {

	using LadderType = Ladder<ComplexOrRealType, InputType>;

public:

	enum
	{
		DIRECTION_X = LadderType::DIRECTION_X,
		DIRECTION_Y = LadderType::DIRECTION_Y,
		DIRECTION_XPY,
		DIRECTION_XMY
	};

	LadderX() { }

	LadderX(SizeType linSize, InputType& io)
	    : ladder_(linSize, io)
	    , linSize_(linSize)
	    , leg_(ladder_.leg())
	{ }

	SizeType maxConnections() const override { return leg_ + 1; }

	SizeType dirs() const override { return 4; }

	SizeType length(SizeType) const override { return this->unimplemented("length"); }

	SizeType translate(SizeType, SizeType, SizeType) const override
	{
		return this->unimplemented("translate");
	}

	SizeType getVectorSize(SizeType dirId) const override
	{
		switch (dirId) {
		case DIRECTION_XPY:
			return linSize_ - leg_;
		case DIRECTION_XMY:
			return linSize_ - leg_;
		}

		return ladder_.getVectorSize(dirId);
	}

	bool connected(SizeType i1, SizeType i2) const override
	{
		if (i1 == i2)
			return false;

		if (ladder_.connected(i1, i2))
			return true;

		SizeType lx          = linSize_ / leg_;
		bool     isPeriodicY = ladder_.isPeriodicY();
		SizeType c1          = i1 / leg_;
		SizeType c2          = i2 / leg_;
		SizeType r1          = i1 % leg_;
		SizeType r2          = i2 % leg_;

		if (c1 == c2)
			return this->neighbors(r1, r2, isPeriodicY, leg_);

		if (r1 == r2)
			return this->neighbors(c1, c2, false, 0);

		return (this->neighbors(r1, r2, isPeriodicY, leg_)
		        && this->neighbors(c1, c2, isPeriodicY, lx));
	}

	// assumes i1 and i2 are connected
	SizeType calcDir(SizeType i1, SizeType i2) const override
	{
		if (ladder_.sameColumn(i1, i2))
			return DIRECTION_Y;
		if (ladder_.sameRow(i1, i2))
			return DIRECTION_X;
		SizeType imin = (i1 < i2) ? i1 : i2;
		if (imin & 1)
			return DIRECTION_XPY;
		return DIRECTION_XMY;
	}

	bool fringe(SizeType i, SizeType smax, SizeType emin) const override
	{
		bool a = (i < emin && i >= smax - 1);
		bool b = (i > smax && i <= emin + 1);
		if (smax & 1)
			return (a || b);
		a = (i < emin && i >= smax - 2);
		b = (i > smax && i <= emin + 2);
		return (a || b);
	}

	// assumes i1 and i2 are connected
	SizeType handle(SizeType i1, SizeType i2) const override
	{
		SizeType dir  = calcDir(i1, i2);
		SizeType imin = (i1 < i2) ? i1 : i2;
		switch (dir) {
		case DIRECTION_X:
			return imin;
		case DIRECTION_Y:
			return imin - imin / leg_;
		case DIRECTION_XPY: // only checked for leg_=2
			return (imin - 1) / leg_;
		case DIRECTION_XMY: // only checked for leg_=2
			return imin / leg_;
		}
		throw RuntimeError("handle: Unknown direction\n");
	}

	// siteNew2 is fringe in the environment
	SizeType getSubstituteSite(SizeType smax, SizeType emin, SizeType siteNew2) const override
	{
		return smax + siteNew2 - emin + 1;
	}

	String label() const override { return "ladderx"; }

	SizeType findReflection(SizeType) const override
	{
		throw RuntimeError("findReflection: unimplemented (sorry)\n");
	}

private:

	LadderType ladder_; // owner
	SizeType   linSize_;
	SizeType   leg_;
}; // class LadderBath
} // namespace PsimagLite

/*@}*/
#endif // GEOMETRY_H
