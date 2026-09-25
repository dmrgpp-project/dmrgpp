// SPDX-FileCopyrightText: Copyright (c) 2009-2014, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file HilbertSpaceFeAs.h
 *
 *  This class represents the Hilbert space for the FeAs Model
 *  States are represented with binary numbers. N bits per site
 *  Bits meaning:
 *  0.....0  empty state
 *  0..1..0  a 1 at location x means state "x"
 *  0..010..010..0, a 1 at location x and a 1 at location y means
 *  2 electrons on the site with states x and y respectively
 *  ...
 *  1111...111111  all ones means N electrons each with a different state
 *
 *  Note: this is a static class
 *
 */
#ifndef HILBERTSPACEFEAS_HEADER_H
#define HILBERTSPACEFEAS_HEADER_H

namespace Dmrg {

//! A class to operate on n-ary numbers (base n)
template <typename Word> class HilbertSpaceFeAs {

	static SizeType orbitals_;

public:

	using HilbertState = Word;

	enum
	{
		SPIN_UP   = 0,
		SPIN_DOWN = 1
	};

	static void setOrbitals(SizeType orbitals) { orbitals_ = orbitals; }

	// Get electronic state on site "j" in binary number "a"
	static Word get(Word const& a, SizeType j)
	{
		SizeType dofs = 2 * orbitals_;
		SizeType k    = dofs * j;
		SizeType ones = (1 << (dofs)) - 1;
		Word     mask = (ones << k);

		mask &= a;
		mask >>= k;
		return mask;
	}

	// Create electron with internal dof "sigma" on site "j" in binary number "a"
	static void create(Word& a, SizeType j, SizeType sigma)
	{
		SizeType dofs = 2 * orbitals_;
		SizeType k    = dofs * j;
		Word     mask = (1 << (k + sigma));
		a |= mask;
	}

	// Destroy electron with internal dof "sigma" on site "j" in binary number "a"
	static void destroy(Word& a, SizeType j, SizeType sigma)
	{
		SizeType dofs = 2 * orbitals_;
		SizeType k    = dofs * j;
		Word     mask = (1 << (k + sigma));
		a &= (~mask);
	}

	// Is there an electron with internal dof "sigma" on site "i" in binary number "ket"?
	static bool isNonZero(Word const& ket, SizeType i, SizeType sigma)
	{

		Word tmp = get(ket, i);
		if (tmp & (1 << sigma))
			return true;

		return false;
	}

	//! returns the number of electrons of internal dof "value" in binary number "data"
	static int getNofDigits(Word const& data, SizeType value)
	{
		SizeType dofs  = 2 * orbitals_;
		int      ret   = 0;
		Word     data2 = data;
		SizeType i     = 0;
		do {
			if ((data & (1 << (dofs * i + value))))
				ret++;
			i++;
		} while (data2 >>= dofs);

		return ret;
	}

	//! Number of electrons with spin spin (sums over bands and sites)
	static int electronsWithGivenSpin(Word const& data, SizeType spin)
	{
		SizeType sum   = 0;
		Word     data2 = data;
		SizeType digit = 0;
		while (data2 > 0) {
			SizeType sigma   = digit % (2 * orbitals_);
			SizeType spin2   = static_cast<SizeType>(sigma / orbitals_);
			Word     thisbit = (data2 & 1);
			data2 >>= 1;
			digit++;
			if (spin == spin2)
				sum += thisbit;
		}

		return sum;
	}

	//! Number of electrons in binary number "data" (sum over all bands)
	static SizeType electrons(const Word& data)
	{
		SizeType sum  = 0;
		SizeType dofs = 2 * orbitals_;
		for (SizeType sector = 0; sector < dofs; sector++)
			sum += calcNofElectrons(data, sector);

		return sum;
	}
	//! Number of electrons with dof sector between i and j excluding
	//! i and j in binary number "ket"
	//!  intended for when i<j
	static int calcNofElectrons(Word const& ket, SizeType i, SizeType j, SizeType sector)
	{
		SizeType dofs = 2 * orbitals_;
		SizeType ii   = i + 1;
		if (ii >= j)
			return 0;
		Word m = 0;
		for (SizeType k = dofs * ii; k < dofs * j; k++)
			m |= (1 << k);
		m = m & ket;
		return getNofDigits(m, sector);
	}

	//! Number of electrons with dof sector in binary number "ket"
	static int calcNofElectrons(Word const& ket, SizeType sector)
	{
		SizeType dofs  = 2 * orbitals_;
		Word     ket2  = ket;
		SizeType digit = 0;
		int      sum   = 0;
		while (ket2 > 0) {
			SizeType sector2 = digit % dofs;
			SizeType thisbit = (ket2 & 1);
			if (sector == sector2)
				sum += thisbit;
			digit++;
			ket2 >>= 1;
		}

		return sum;
	}

}; // class HilbertSpaceFeAs

template <typename Word> SizeType HilbertSpaceFeAs<Word>::orbitals_ = 2;
} // namespace Dmrg

/*@}*/
#endif
