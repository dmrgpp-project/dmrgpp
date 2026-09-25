// SPDX-FileCopyrightText: Copyright (c) 2009-2014, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file HilbertSpaceImmm.h
 *
 *  This class represents the Hilbert space for the Immm Model
 *  States are represented with binary numbers. N(site) bits per site
 *  Bits meaning:
 *  0.....0  empty state
 *  0..1..0  a 1 at location x means state "x"
 *  0..010..010..0, a 1 at location x and a 1 at location y means 2
 *  electrons on the site with states x and y respectively
 *  ...
 *  1111...111111  all ones means N electrons each with a different state
 *
 *  Note: this is a class
 *  Note: Length of state depends on site
 *
 */
#ifndef HILBERTSPACE_IMMM_H
#define HILBERTSPACE_IMMM_H
#include <iostream>
#include <vector>

namespace Dmrg {

//! A class to operate on n-ary numbers (base n)
template <typename Word> class HilbertSpaceImmm {

public:

	using HilbertState = Word;

	static const SizeType NUMBER_OF_SPINS = 2;
	enum
	{
		SPIN_UP   = 0,
		SPIN_DOWN = 1
	};

	HilbertSpaceImmm(SizeType maxOrbitals)
	    : maxOrbitals_(maxOrbitals)
	{ }

	template <typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType& mres, SizeType, PsimagLite::String msg = "") const
	{
		PsimagLite::String str = msg;
		str += "HilbertSpaceImmm";

		mres.memResolv(&maxOrbitals_, sizeof(*this), str + " maxOrbitals");

		return sizeof(*this);
	}

	SizeType dOf() const { return 2 * maxOrbitals_; }

	// Get electronic state on site "j" in binary number "a"
	Word get(Word const& a, SizeType j) const
	{
		SizeType k    = degreesOfFreedomUpTo(j);
		SizeType ones = (1 << (dOf())) - 1;
		Word     mask = (ones << k);

		mask &= a;
		mask >>= k;
		return mask;
	}

	// Create electron with internal dof  "sigma" on site "j" in binary number "a"
	void create(Word& a, SizeType j, SizeType sigma) const
	{
		SizeType k    = degreesOfFreedomUpTo(j);
		Word     mask = (1 << (k + sigma));
		a |= mask;
	}

	// Is there an electron with internal dof  "sigma" on site "i" in binary number "ket"?
	bool isNonZero(Word const& ket, SizeType i, SizeType sigma) const
	{

		Word tmp = get(ket, i);
		if (tmp & (1 << sigma))
			return true;

		return false;
	}

	//! returns the number of electrons of internal dof "value" in binary number "data"
	int getNofDigits(const Word& data, SizeType value) const
	{
		int      ret   = 0;
		Word     data2 = data;
		SizeType i     = 0;
		SizeType dof   = 0;

		do {
			SizeType k = degreesOfFreedomUpTo(i);
			dof        = dOf();
			if ((data & (1 << (k + value))))
				ret++;
			i++;
		} while (data2 >>= dof);

		return ret;
	}

	//! Number of electrons with spin spin (sums over bands)
	int electronsWithGivenSpin(Word const& data, SizeType, SizeType spin) const
	{

		SizeType norb   = dOf() / NUMBER_OF_SPINS;
		SizeType beginX = spin * norb;
		SizeType endX   = beginX + norb;
		SizeType sum    = 0;

		for (SizeType x = beginX; x < endX; x++)
			sum += getNofDigits(data, x);

		return sum;
	}

	//! Number of electrons for data
	SizeType electrons(const Word& data) const
	{
		SizeType sum = 0;
		for (SizeType sector = 0; sector < dOf(); sector++)
			sum += calcNofElectrons(data, 0, sector);

		return sum;
	}

	//! Number of electrons with dof sector between i and j
	//! excluding i and j in binary number "ket"
	//!  intended for when i<j
	int calcNofElectrons(const Word& ket, SizeType i, SizeType j, SizeType sector) const
	{
		SizeType ii = i + 1;
		if (ii >= j)
			return 0;
		Word m = 0;
		for (SizeType site = ii; site < j; site++) {
			SizeType k   = degreesOfFreedomUpTo(site);
			SizeType dof = dOf();
			for (SizeType sigma = 0; sigma < dof; sigma++)
				m |= (1 << (k + sigma));
		}
		m &= ket;
		return getNofDigits(m, sector);
	}

	//! Number of electrons with dof sector on site i in binary number "ket"
	int calcNofElectrons(Word const& ket, SizeType i, SizeType sector) const
	{
		Word m = 0;

		SizeType k   = degreesOfFreedomUpTo(i);
		SizeType dof = dOf();
		for (SizeType sigma = 0; sigma < dof; sigma++)
			m |= (1 << (k + sigma));

		m &= ket;
		return getNofDigits(m, sector);
	}

	void write(PsimagLite::String label, PsimagLite::IoNg::Out::Serializer& io) const
	{
		PsimagLite::String label2 = label + "/HilbertSpaceImmm";
		io.createGroup(label2);
		io.write(label2 + "/maxOrbitals_", maxOrbitals_);
	}

private:

	SizeType degreesOfFreedomUpTo(SizeType j) const { return dOf() * j; }

	// serializr start class HilbertSpaceImmm
	// serializr normal maxOrbitals_
	SizeType maxOrbitals_;
}; // class HilbertSpaceImmm
} // namespace Dmrg

/*@}*/
#endif // HILBERTSPACE_IMMM_H
