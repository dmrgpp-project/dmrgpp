// SPDX-FileCopyrightText: Copyright (c) 2009, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 2.0.0]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file SpinSquared.h
 *
 *  Encapsulates the calculation of S^2 from creation
 *  and anihilation operators. Works in conjunction with SpinSquaredHelper
 *
 */
#ifndef SPIN_SQUARED_H
#define SPIN_SQUARED_H

namespace Dmrg {
template <typename CallbackType> class SpinSquared {

	using Word      = typename CallbackType::Word;
	using FieldType = typename CallbackType::FieldType;

public:

	enum
	{
		SPIN_UP   = 0,
		SPIN_DOWN = 1
	};

	enum
	{
		ORBITAL_A = 0,
		ORBITAL_B = 1
	};

	SpinSquared(CallbackType& callback, int NUMBER_OF_ORBITALS1, int DEGREES_OF_FREEDOM1)
	    : callback_(callback)
	    , NUMBER_OF_ORBITALS(NUMBER_OF_ORBITALS1)
	    , DEGREES_OF_FREEDOM(DEGREES_OF_FREEDOM1)
	{ }

	template <typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType& mres, SizeType, PsimagLite::String msg = "") const
	{
		PsimagLite::String str = msg;
		str += "SpinSquared";

		const char* start = reinterpret_cast<const char*>(this);
		const char* end   = reinterpret_cast<const char*>(&NUMBER_OF_ORBITALS);
		SizeType    total = end - start;
		mres.push(SomeMemResolvType::MEMORY_HEAPPTR,
		          SomeMemResolvType::SIZEOF_HEAPREF,
		          this,
		          str + " ref to callback");

		mres.memResolv(&callback_, 0, str + " callback");

		start = end;
		end   = reinterpret_cast<const char*>(&DEGREES_OF_FREEDOM);
		total += mres.memResolv(
		    &NUMBER_OF_ORBITALS, end - start, str + " NUMBER_OF_ORBITALS");

		total += mres.memResolv(
		    &DEGREES_OF_FREEDOM, sizeof(*this) - total, str + " DEGREES_OF_FREEDOM");

		return total;
	}

	void doOnePairOfSitesA(const Word& ket, SizeType i, SizeType j) const
	{
		FieldType value = 0.5;
		Word      bra   = ket;
		int       ret   = c(bra, j, ORBITAL_A, SPIN_UP);
		ret             = cDagger(bra, j, ORBITAL_A, SPIN_DOWN, ret);
		ret             = c(bra, i, ORBITAL_A, SPIN_DOWN, ret);
		ret             = cDagger(bra, i, ORBITAL_A, SPIN_UP, ret);
		if (ret == 0)
			callback_(ket, bra, value);
		if (NUMBER_OF_ORBITALS == 1)
			return;

		bra = ket;
		ret = c(bra, j, ORBITAL_B, SPIN_UP);
		ret = cDagger(bra, j, ORBITAL_B, SPIN_DOWN, ret);
		ret = c(bra, i, ORBITAL_A, SPIN_DOWN, ret);
		ret = cDagger(bra, i, ORBITAL_A, SPIN_UP, ret);
		if (ret == 0)
			callback_(ket, bra, value);

		bra = ket;
		ret = c(bra, j, ORBITAL_A, SPIN_UP);
		ret = cDagger(bra, j, ORBITAL_A, SPIN_DOWN, ret);
		ret = c(bra, i, ORBITAL_B, SPIN_DOWN, ret);
		ret = cDagger(bra, i, ORBITAL_B, SPIN_UP, ret);
		if (ret == 0)
			callback_(ket, bra, value);

		bra = ket;
		ret = c(bra, j, ORBITAL_B, SPIN_UP);
		ret = cDagger(bra, j, ORBITAL_B, SPIN_DOWN, ret);
		ret = c(bra, i, ORBITAL_B, SPIN_DOWN, ret);
		ret = cDagger(bra, i, ORBITAL_B, SPIN_UP, ret);
		if (ret == 0)
			callback_(ket, bra, value);
	}

	void doOnePairOfSitesB(const Word& ket, SizeType i, SizeType j) const
	{
		FieldType value = 0.5;
		Word      bra   = ket;
		int       ret   = c(bra, j, ORBITAL_A, SPIN_DOWN);
		ret             = cDagger(bra, j, ORBITAL_A, SPIN_UP, ret);
		ret             = c(bra, i, ORBITAL_A, SPIN_UP, ret);
		ret             = cDagger(bra, i, ORBITAL_A, SPIN_DOWN, ret);
		if (ret == 0)
			callback_(ket, bra, value);
		if (NUMBER_OF_ORBITALS == 1)
			return;

		bra = ket;
		ret = c(bra, j, ORBITAL_B, SPIN_DOWN);
		ret = cDagger(bra, j, ORBITAL_B, SPIN_UP, ret);
		ret = c(bra, i, ORBITAL_A, SPIN_UP, ret);
		ret = cDagger(bra, i, ORBITAL_A, SPIN_DOWN, ret);
		if (ret == 0)
			callback_(ket, bra, value);

		bra = ket;
		ret = c(bra, j, ORBITAL_A, SPIN_DOWN);
		ret = cDagger(bra, j, ORBITAL_A, SPIN_UP, ret);
		ret = c(bra, i, ORBITAL_B, SPIN_UP, ret);
		ret = cDagger(bra, i, ORBITAL_B, SPIN_DOWN, ret);
		if (ret == 0)
			callback_(ket, bra, value);

		bra = ket;
		ret = c(bra, j, ORBITAL_B, SPIN_DOWN);
		ret = cDagger(bra, j, ORBITAL_B, SPIN_UP, ret);
		ret = c(bra, i, ORBITAL_B, SPIN_UP, ret);
		ret = cDagger(bra, i, ORBITAL_B, SPIN_DOWN, ret);
		if (ret == 0)
			callback_(ket, bra, value);
	}

	void doDiagonal(const Word& ket, SizeType i, SizeType j) const
	{
		FieldType value = spinZ(ket, i) * spinZ(ket, j);
		callback_(ket, ket, value);
	}

	FieldType spinZ(const Word& ket, SizeType i) const
	{
		int sum = 0;
		for (SizeType gamma = 0; gamma < SizeType(NUMBER_OF_ORBITALS); gamma++) {
			sum += n(ket, i, gamma, SPIN_UP);
			sum -= n(ket, i, gamma, SPIN_DOWN);
		}
		return sum * 0.5;
	}

	void write(PsimagLite::String, PsimagLite::IoNg::Out::Serializer&) const { }

private:

	int c(Word& ket, SizeType i, SizeType gamma, SizeType spin, int ret = 0) const
	{
		if (ret < 0)
			return ret;
		SizeType shift_ = (DEGREES_OF_FREEDOM * i);
		Word     mask   = 1 << (shift_ + gamma + spin * NUMBER_OF_ORBITALS);
		if ((ket & mask) == 0)
			return -1;
		ket ^= mask;
		return 0;
	}

	int cDagger(Word& ket, SizeType i, SizeType gamma, SizeType spin, int ret = 0) const
	{
		if (ret < 0)
			return ret;
		SizeType shift_ = (DEGREES_OF_FREEDOM * i);
		Word     mask   = 1 << (shift_ + gamma + spin * NUMBER_OF_ORBITALS);
		if ((ket & mask) > 0)
			return -1;
		ket |= mask;
		return 0;
	}

	int n(const Word& ket, SizeType i, SizeType gamma, SizeType spin) const
	{
		SizeType shift_ = (DEGREES_OF_FREEDOM * i);
		Word     mask   = 1 << (shift_ + gamma + spin * NUMBER_OF_ORBITALS);
		if ((ket & mask) > 0)
			return 1;
		return 0;
	}

	// serializr start class SpinSquared
	// serializr ref callback_ this
	CallbackType& callback_;
	// serializr normal NUMBER_OF_ORBITALS
	int const NUMBER_OF_ORBITALS;
	// serializr normal DEGREES_OF_FREEDOM
	int const DEGREES_OF_FREEDOM;
}; // SpinSquared
} // namespace Dmrg

/*@}*/
#endif
