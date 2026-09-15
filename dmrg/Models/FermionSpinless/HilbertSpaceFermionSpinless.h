// SPDX-FileCopyrightText: Copyright (c) 2009-2015, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file HilbertSpaceFermionSpinless.h
 *
 *  This class represents the Hilbert space for the FermionSpinless Model
 *  States are represented with binary numbers. One bit per site
 *
 *  Note: this is a static class
 *
 */
#ifndef DMRG_HILBERTSPACE_FERMIONSPINLESS_H
#define DMRG_HILBERTSPACE_FERMIONSPINLESS_H

namespace Dmrg {

template <typename Word> class HilbertSpaceFermionSpinless {
public:

	using HilbertState = Word;

	//! For state "a" set electron on site "j" to value "value"
	static void set(Word& a, int j, int value)
	{
		Word mask;
		switch (value) {
		case 0:
			mask = (1 << j);
			a &= (~mask);
			return;
		case 1:
			mask = (1 << j);
			a |= mask;
			return;
		default:
			std::cerr << "value=" << value << "\n";
			throw PsimagLite::RuntimeError("set: invalid value.\n");
		}
	}

	// Get electronic state on site "j" in binary number "a"
	static int get(Word const& a, int j)
	{
		Word mask = (1 << j);
		mask &= a;
		mask >>= j;
		assert(mask <= 1);
		return mask;
	}

	// Destroy electron with internal dof  "sigma" on site "j" in binary number "a"
	static void destroy(Word& a, int j, int sigma)
	{
		assert(sigma == 0);
		Word mask;
		switch (sigma) {
		case 0:
			mask = (1 << j);
			a &= (~mask);
			return;
		default:
			std::cerr << "sigma=" << sigma << "\n";
			throw PsimagLite::RuntimeError("destroy: invalid value.\n");
		}
	}

	// Create electron with internal dof  "sigma" on site "j" in binary number "a"
	static void create(Word& a, int j, int sigma)
	{
		assert(sigma == 0);
		Word mask;
		switch (sigma) {
		case 0:
			mask = (1 << j);
			a |= mask;
			return;
		default:
			std::cerr << "sigma=" << sigma << "\n";
			throw PsimagLite::RuntimeError("create: invalid value.\n");
		}
	}

	// Is there an electron with internal dof
	// "sigma" on site "i" in binary number "ket"?
	static bool isNonZero(const Word& ket, int i, int sigma)
	{
		assert(sigma == 0);
		int tmp = get(ket, i);
		if ((tmp & 1) && sigma == 0)
			return true;

		return false;
	}

	// returns the number of electrons of internal dof "value" in binary number "data"
	static int getNofDigits(Word const& data, int value)
	{
		assert(value == 0);
		int  ret   = 0;
		Word data2 = data;
		int  i     = 0;
		do {
			if ((data & (1 << (i + value))))
				ret++;
			i++;
		} while (data2 >>= 1);

		return ret;
	}

	// Number of electrons with dof sector between i and
	// j excluding i and j in binary number "ket"
	//  intended for when i<j
	static int calcNofElectrons(Word const& ket, int i, int j, int sector)
	{
		int ii = i + 1;
		if (ii >= j)
			return 0;
		Word     m   = 0;
		SizeType end = j;
		for (SizeType k = ii; k < end; k++)
			m |= (1 << k);
		m = m & ket;
		return getNofDigits(m, sector);
	}

}; // class HilbertSpaceFermionSpinless
} // namespace Dmrg

/*@}*/
#endif
