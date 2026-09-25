// SPDX-FileCopyrightText: Copyright (c) 2009-2014, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file HilbertSpaceHubbardHolstein.h
 *
 *  This class represents the Hilbert space for the Hubbard Holstein Model
 *  States are represented with binary numbers. N bits per site
 * (2 for the fermions, N-2 for the phonons)
 *  Bits meaning:
 *  0....00  empty state
 *  0..1.00  a 1 at location x>2 means state with N-2 maxphonons and zero fermions
 *  1111111  all ones means 2 electrons
 *
 *  Note: this is a static class
 *
 */
#ifndef HILBERTSPACE_HUBBARD_HOLSTEIN_HEADER_H
#define HILBERTSPACE_HUBBARD_HOLSTEIN_HEADER_H

#include "Utils.h"

namespace Dmrg {

//! A class to operate on n-ary numbers (base n)
template <typename Word> class HilbertSpaceHubbardHolstein {

public:

	using HilbertState = Word;

	enum class SpinEnum
	{
		SPIN_UP,
		SPIN_DOWN
	};

	// Get electronic state from combined electron and phonon ket a
	static Word getF(Word a) { return (a & 3); }

	// Get phononic state from combined electron and phonon ket a
	static Word getP(Word a) { return (a >> 2); }

	// Create electron with internal dof "sigma" in binary number "a"
	static void createF(Word& a, SizeType sigma)
	{
		Word mask = (1 << sigma);
		assert((a & mask) == 0);
		a |= mask;
	}

	// Create phonon
	static void createP(Word& a)
	{
		const SizeType nphonons = getP(a);
		const Word     stateP   = 1 + nphonons;
		const Word     maskP    = (stateP << 2);
		const Word     maskF    = getF(a);
		assert(maskF < 4);
		a = (maskP | maskF);
	}

	static SizeType electronsWithGivenSpin(Word a, SpinEnum spin)
	{
		const SizeType f = getF(a);
		return (spin == SpinEnum::SPIN_UP) ? (f & 1) : (f >> 1);
	}
}; // class HilbertSpaceHubbardHolstein

} // namespace Dmrg
/*@}*/
#endif
