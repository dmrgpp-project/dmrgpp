// SPDX-FileCopyrightText: Copyright (c) 2009-2014, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file HilbertSpaceHubbardHolstein.h
 *
 *  This class represents the Hilbert space for the Hubbard Holstein Model Spinless
 *  States are represented with binary numbers. N bits per site
 * (1 for the fermions, N-1 for the phonons)
 *  Bits meaning:
 *  0....00  empty state
 *  0..1.00  a 1 at location x>1 means state with N-1 maxphonons and zero fermions
 *  1111111  all ones means 1 electrons
 *
 *  Note: this is a static class
 *
 */
#ifndef HILBERTSPACE_HUBBARD_HOLSTEIN_SPINLESS_HEADER_H
#define HILBERTSPACE_HUBBARD_HOLSTEIN_SPINLESS_HEADER_H

#include "Utils.h"

namespace Dmrg {

//! A class to operate on n-ary numbers (base n)
template <typename Word> class HilbertSpaceHubbardHolsteinSpinless {

public:

	using HilbertState = Word;

	// Get electronic state from combined electron and phonon ket a
	static Word getF(Word a) { return (a & 1); }

	// Get phononic state from combined electron and phonon ket a
	static Word getP(Word a) { return (a >> 1); }

	// Create electron in binary number "a"
	static void createF(Word& a)
	{
		Word mask = 1;
		assert((a & mask) == 0);
		a |= mask;
	}

	// Create phonon
	static void createP(Word& a)
	{
		const SizeType nphonons = getP(a);
		const Word     stateP   = 1 + nphonons;
		const Word     maskP    = (stateP << 1);
		const Word     maskF    = getF(a);
		assert(maskF < 2);
		a = (maskP | maskF);
	}

	static SizeType electronsWithGivenSpin(Word a)
	{
		const SizeType f = getF(a);
		return (f & 1);
	}
}; // class HilbertSpaceHubbardHolstein

} // namespace Dmrg
/*@}*/
#endif
