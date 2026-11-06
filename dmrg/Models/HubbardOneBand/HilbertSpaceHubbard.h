// SPDX-FileCopyrightText: Copyright (c) 2009-2013, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 2.0.0]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file HilbertSpaceHubbard.h
 *
 *  This class represents the Hilbert space for the Hubbard Model
 *  States are represented with binary numbers. Two bits per site
 *  Bits meaning:
 *  00 0 empty state
 *  01 1 up electron
 *  10 2 down electron
 *  11 3 doubly-occupied state
 *
 *  Note: this is a static class
 *
 */
#ifndef HILBERTSPACEHUBBARD_HEADER_H
#define HILBERTSPACEHUBBARD_HEADER_H

namespace Dmrg {

//! A class to operate on quaterny numbers (base 4)
template <typename Word> class HilbertSpaceHubbard {
public:

	using HilbertState         = Word;
	static int const SPIN_UP   = 0;
	static int const SPIN_DOWN = 1;

	//! For state "a" set electron on site "j" to value "value"
	static void set(Word& a, int j, int value)
	{
		Word mask;
		switch (value) {
		case 0:
			mask = (1 << (2 * j)) | (1 << (2 * j + 1));
			a &= (~mask);
			return;
		case 1:
			mask = (1 << (2 * j));
			a |= mask;
			mask = (1 << (2 * j + 1));
			a &= (~mask);
			return;
		case 2:
			mask = (1 << (2 * j + 1));
			a |= mask;
			mask = (1 << (2 * j));
			a &= (~mask);
			return;
		case 3:
			mask = (1 << (2 * j)) | (1 << (2 * j + 1));
			a |= (mask);
			return;
		default:
			std::cerr << "value=" << value << "\n";
			err("HilbertSpaceHubbard: set: invalid value.\n");
		}
	}

	// Get electronic state on site "j" in binary number "a"
	static int get(Word const& a, int j)
	{
		Word mask = (1 << (2 * j)) | (1 << (2 * j + 1));
		mask &= a;
		mask >>= (2 * j);
		if (mask > 3) {
			std::cerr << "Error: mask=" << mask << "\n";
			err("HilbertSpaceHubbard: get: invalid  mask.\n");
		}

		return mask;
	}

	// Destroy electron with internal dof  "sigma" on site "j" in binary number "a"
	static void destroy(Word& a, int j, int sigma)
	{
		Word mask;
		switch (sigma) {
		case 0:
			mask = (1 << (2 * j));
			a &= (~mask);
			return;
		case 1:
			mask = (1 << (2 * j + 1));
			a &= (~mask);
			return;
		default:
			std::cerr << "sigma=" << sigma << "\n";
			err("HilbertSpaceHubbard: destroy: invalid value.\n");
		}
	}

	// Create electron with internal dof  "sigma" on site "j" in binary number "a"
	static void create(Word& a, int j, int sigma)
	{
		Word mask;
		switch (sigma) {
		case 0:
			mask = (1 << (2 * j));
			a |= mask;
			return;
		case 1:
			mask = (1 << (2 * j + 1));
			a |= mask;
			return;
		default:
			std::cerr << "sigma=" << sigma << "\n";
			err("HilbertSpaceHubbard: create: invalid value.\n");
		}
	}

	// Is there an electron with internal dof  "sigma" on site "i" in binary number "ket"?
	static bool isNonZero(Word const& ket, int i, int sigma)
	{
		int          tmp  = get(ket, i);
		HilbertState mask = (1 << sigma);
		return (tmp & mask);
	}

	//! returns the number of electrons of internal dof "value" in binary number "data"
	static int getNofDigits(Word const& data, int value)
	{
		int  ret   = 0;
		Word data2 = data;
		int  i     = 0;
		do {
			if ((data & (1 << (2 * i + value))))
				ret++;
			i++;
		} while (data2 >>= 1);

		return ret;
	}

	// Number of electrons with dof sector between i and j
	// excluding i and j in binary number "ket"
	//  intended for when i<j
	static int calcNofElectrons(Word const& ket, int i, int j, int sector)
	{
		int ii = i + 1;
		if (ii >= j)
			return 0;
		Word     m   = 0;
		SizeType end = 2 * j;
		for (SizeType k = 2 * ii; k < end; k++)
			m |= (1 << k);
		m = m & ket;
		return getNofDigits(m, sector);
	}
}; // class HilbertSpaceHubbard
} // namespace Dmrg

/*@}*/
#endif
