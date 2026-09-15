// SPDX-FileCopyrightText: Copyright (c) 2009-2013-2018, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 2.]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup PsimagLite */
/*@{*/

/*! \file Honeycomb.h
 *
 *  HoneyComb lattice implementation (started March 2018)
 */
#ifndef PSI_GEOMETRY_HONEYCOMB_H
#define PSI_GEOMETRY_HONEYCOMB_H

#include "GeometryBase.h"
#include <cstdlib>
#include <stdexcept>

namespace PsimagLite {

template <typename ComplexOrRealType, typename InputType>
class Honeycomb : public GeometryBase<ComplexOrRealType, InputType> {

public:

	enum Dir
	{
		DIR_X = 0,
		DIR_Y = 1,
		DIR_Z = 2
	};
	enum GeLe
	{
		GREATER_OR_EQUAL,
		LESS_OR_EQUAL
	};

	Honeycomb(SizeType linSize, InputType& io)
	    : linSize_(linSize)
	    , periodicY_(false)
	    , periodicX_(false)
	{
		io.readline(ly_, "HoneycombLy=");
		int x = 0;
		io.readline(x, "IsPeriodicY=");
		periodicY_ = (x > 0);

		try {
			io.readline(x, "IsPeriodicX=");
			periodicX_ = (x > 0);
		} catch (std::exception&) { }
	}

	SizeType maxConnections() const override
	{
		// max connections broken when splitting the lattice
		return ly_ + 1; // overestimate probably
	}

	SizeType dirs() const override { return 3; }

	SizeType length(SizeType) const override { return this->unimplemented("length"); }

	SizeType translate(SizeType, SizeType, SizeType) const override
	{
		return this->unimplemented("translate");
	}

	SizeType getVectorSize(SizeType dirId) const override
	{
		throw RuntimeError("Honeycomb::getVectorSize() unimplemented\n");
	}

	bool connected(SizeType i1, SizeType i2) const override
	{
		return connectedInternal(i1, i2).first;
	}

	// assumes i1 and i2 are connected
	SizeType calcDir(SizeType i1, SizeType i2) const override
	{
		std::pair<bool, Dir> bdir = connectedInternal(i1, i2);
		assert(bdir.first);
		return bdir.second;
	}

	bool fringe(SizeType i, SizeType smax, SizeType emin) const override
	{
		if (i <= smax)
			return isThereAneighbor(i, emin, linSize_);

		if (i >= emin)
			return isThereAneighbor(i, 0, smax);

		return false;
	}

	// assumes i1 and i2 are connected
	SizeType handle(SizeType i1, SizeType i2) const override
	{
		throw RuntimeError("Honeycomb::handle() unimplemented\n");
	}

	// siteNew2 is fringe in the environment
	SizeType getSubstituteSite(SizeType smax, SizeType emin, SizeType siteNew2) const override
	{
		throw RuntimeError("Honeycomb::getSubstituteSite() unimplemented\n");
	}

	String label() const override { return "Honeycomb"; }

	SizeType findReflection(SizeType) const override
	{
		throw RuntimeError("findReflection: unimplemented (sorry)\n");
	}

	template <class Archive> void write(Archive&, const unsigned int) { }

private:

	std::pair<bool, Dir> connectedInternal(SizeType ii1, SizeType ii2) const
	{
		std::pair<bool, Dir> falseDir(false, DIR_X);
		if (ii1 == ii2)
			return falseDir;

		bool     normal = (ii1 < ii2);
		SizeType i1     = (normal) ? ii1 : ii2;
		SizeType i2     = (normal) ? ii2 : ii1;

		SizeType x1 = 0;
		SizeType y1 = 0;
		getCoordinates(x1, y1, i1);

		SizeType x2 = 0;
		SizeType y2 = 0;
		getCoordinates(x2, y2, i2);

		if (isDirectionX(x1, y1, x2, y2))
			return std::pair<bool, Dir>(true, DIR_X);

		if (isDirectionY(x1, y1, x2, y2))
			return std::pair<bool, Dir>(true, DIR_Y);

		if (isDirectionZ(x1, y1, x2, y2))
			return std::pair<bool, Dir>(true, DIR_Z);

		return falseDir;
	}

	void getCoordinates(SizeType& x, SizeType& y, SizeType i) const
	{
		div_t q = std::div(static_cast<int>(i), static_cast<int>(ly_));
		x       = q.quot;
		y       = q.rem;
	}

	// assumes i1 < i2
	bool isDirectionZ(SizeType x1, SizeType y1, SizeType x2, SizeType y2) const
	{
		if (x1 != x2)
			return false;
		assert(y1 < y2);
		SizeType d = y2 - y1;
		if ((d == 1) && (y1 & 1))
			return true;

		if (!periodicY_)
			return false;

		// periodic in y direction
		return (d + 1 == ly_);
	}

	// assumes i1 < i2
	bool isDirectionY(SizeType x1, SizeType y1, SizeType x2, SizeType y2) const
	{
		if (x1 == x2)
			return isDirectionYsameCol(y1, y2);

		SizeType dy = y2 - y1;
		if (dy != 1)
			return false;

		assert(x1 < x2);
		SizeType dx = x2 - x1;
		if ((dx == 1) && isMultipleOf3(y1) && (y2 < y1))
			return true;

		if (!periodicX_)
			return false;

		// periodic in x direction
		SizeType lx = linSize_ / ly_;
		if (dx + 1 != lx)
			return false;

		return (isMultipleOf3(y2) && (y1 < y2));
	}

	// assumes i1 < i2
	bool isDirectionX(SizeType x1, SizeType y1, SizeType x2, SizeType y2) const
	{
		if (x1 == x2)
			return isDirectionXsameCol(y1, y2);

		SizeType dy = y2 - y1;
		if (dy != 1)
			return false;

		assert(x1 < x2);
		SizeType dx = x2 - x1;
		if ((dx == 1) && isMultipleOf4(y1) && (y1 < y2))
			return true;

		if (!periodicX_)
			return false;

		// periodic in x direction
		SizeType lx = linSize_ / ly_;
		if (dx + 1 != lx)
			return false;

		return (isMultipleOf4(y2) && (y2 < y1));
	}

	bool isDirectionYsameCol(SizeType y1, SizeType y2) const
	{
		assert(y1 < y2);
		SizeType dy = y2 - y1;
		if (dy == 1 && isMultipleOf4(y1))
			return true;
		return (periodicY_ && (dy + 1 == ly_) && isMultipleOf4(y2));
	}

	bool isDirectionXsameCol(SizeType y1, SizeType y2) const
	{
		assert(y1 < y2);
		SizeType dy = y2 - y1;
		if (dy == 1 && is4nPlus2(y1))
			return true;
		return (periodicY_ && (dy + 1 == ly_) && is4nPlus2(y2));
	}

	bool isThereAneighbor(SizeType ind, SizeType start, SizeType end) const
	{
		for (SizeType i = start; i < end; ++i) {
			if (connectedInternal(ind, i).first)
				return true;
		}

		return false;
	}

	static bool isMultipleOf3(SizeType y) { return ((y % 3) == 0); }

	static bool isMultipleOf4(SizeType y) { return ((y % 4) == 0); }

	static bool is4nPlus2(SizeType y)
	{
		if (y < 2)
			return false;
		y -= 2;
		return ((y % 4) == 0);
	}

	SizeType linSize_;
	SizeType ly_;
	bool     periodicY_;
	bool     periodicX_;
};
} // namespace PsimagLite

/*@}*/
#endif // GEOMETRY_H
