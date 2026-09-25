// SPDX-FileCopyrightText: Copyright (c) 2009-2014, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file HilbertSpaceIsingMultiOrb.h
 *
 *  This class represents the Hilbert space for the Ising Model
 *  States are represented with binary numbers. N bits per site
 *  Bits meaning:
 *  0.....0  all down spin state
 *  0..1..0  a 1 at location x means spin up on location x
 *  1111...111111  all ones means all up spins
 *
 *  Note: this is a static class
 *
 */
#ifndef HILBERTSPACEISING_HEADER_H
#define HILBERTSPACEISING_HEADER_H

namespace Dmrg {

//! A class to operate on n-ary numbers (base n)
template <typename Word> class HilbertSpaceIsingMultiOrb {

	static SizeType orbitals_;

public:

	using HilbertState = Word;
	using WordType     = unsigned int long;

	enum
	{
		SPIN_UP   = 0,
		SPIN_DOWN = 1
	};

	static void setOrbitals(SizeType orbitals) { orbitals_ = orbitals; }

	static bool isBitZeroAt(SizeType pos, SizeType ket)
	{
		SizeType mask = 1ul;
		bool     flag = true;
		for (SizeType i = 1; i < pos; i++)
			mask <<= 1;

		SizeType check = ket & mask;
		if (check > 0)
			flag = false;

		return flag;
	}

private:

	static PsimagLite::Vector<WordType>::Type bitmask_;

}; // class HilbertSpaceIsingMultiOrb

template <typename Word> SizeType HilbertSpaceIsingMultiOrb<Word>::orbitals_ = 1;
} // namespace Dmrg

/*@}*/
#endif
