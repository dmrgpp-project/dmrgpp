// SPDX-FileCopyrightText: Copyright (c) 2009-2014, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup PsimagLite */
/*@{*/

/*! \file LongChain.h
 *
 *  DOC NEEDED FIXME
 */
#ifndef PSI_GEOM_LONG_CHAIN_H
#define PSI_GEOM_LONG_CHAIN_H
#include "GeometryBase.h"

namespace PsimagLite {

template <typename ComplexOrRealType, typename InputType>
class LongChain : public GeometryBase<ComplexOrRealType, InputType> {

public:

	enum
	{
		DIRECTION_X
	};

	LongChain()
	    : linSize_(0)
	    , isPeriodic_(false)
	    , distance_(1)
	{ }

	LongChain(SizeType linSize, InputType& io)
	    : linSize_(linSize)
	    , isPeriodic_(false)
	    , distance_(1)
	{
		try {
			int x = 0;
			io.readline(x, "IsPeriodicX=");
			isPeriodic_ = (x > 0) ? true : false;
			if (isPeriodic_)
				std::cerr << "LongChain::ctor(): periodic\n";
		} catch (std::exception& e) { }

		try {
			int x = 0;
			io.readline(x, "LongChainDistance=");
			distance_ = x;
		} catch (std::exception& e) { }

		if (linSize_ <= distance_)
			RuntimeError("LongChain::ctor()\n");
	}

	SizeType maxConnections() const override { return (isPeriodic_) ? linSize_ : 1; }

	SizeType dirs() const override { return 1; }

	SizeType handle(SizeType i, SizeType j) const override { return (i < j) ? i : j; }

	SizeType getVectorSize(SizeType dirId) const override
	{
		assert(dirId == DIRECTION_X);
		const SizeType oneOrZero = (isPeriodic_) ? 1 : 0;
		return linSize_ - distance_ + oneOrZero;
	}

	bool connected(SizeType i1, SizeType i2) const override
	{
		if (i1 == i2)
			return false;
		SizeType imin = (i1 < i2) ? i1 : i2;
		SizeType imax = (i1 > i2) ? i1 : i2;
		bool     b    = (imax - imin == distance_);
		if (!isPeriodic_)
			return b;
		bool b2 = (imax - imin == linSize_ - distance_);
		return (b || b2);
	}

	// assumes i1 and i2 are connected
	SizeType calcDir(SizeType, SizeType) const override { return DIRECTION_X; }

	bool fringe(SizeType i, SizeType smax, SizeType emin) const override
	{
		SizeType emin2 = smax + 1;
		if (i <= smax) {
			SizeType iPlus = i + distance_;
			return (iPlus < linSize_ && iPlus >= emin2);
		}

		return (i >= distance_ && (i - distance_) < emin);
	}

	// siteNew2 is fringe in the environment
	SizeType getSubstituteSite(SizeType smax, SizeType emin, SizeType siteNew2) const override
	{
		assert(siteNew2 >= emin);
		SizeType tmp = siteNew2 - emin + smax + 1;
		assert(tmp < linSize_);
		return tmp;
	}

	String label() const override { return "longchain"; }

	SizeType findReflection(SizeType site) const override { return linSize_ - site - 1; }

	SizeType length(SizeType i) const override
	{
		assert(i == 0);
		return linSize_;
	}

	SizeType translate(SizeType site, SizeType dir, SizeType amount) const override
	{
		assert(dir == 0);

		site += amount;
		while (site >= linSize_)
			site -= linSize_;
		return site;
	}

	template <class Archive> void write(Archive&, const unsigned int)
	{
		throw RuntimeError("LongChain::write(): unimplemented\n");
	}

private:

	SizeType linSize_;
	bool     isPeriodic_;
	SizeType distance_;
}; // class LongChain
} // namespace PsimagLite

/*@}*/
#endif // LADDER_H
