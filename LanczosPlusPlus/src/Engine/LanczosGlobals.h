// SPDX-FileCopyrightText: Copyright (c) 2009-2012, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [LanczosPlusPlus, Version 1.0.0]
// [by G.A., Oak Ridge National Laboratory]


/** \ingroup LanczosPlusPlus */
/*@{*/

/*! \file LanczosGlobals.h
 *
 *
 *
 */
#ifndef LANCZOS_PROGRAM_LIMITS_H
#define LANCZOS_PROGRAM_LIMITS_H
#include <PsimagLite/BitManip.h>
#include <PsimagLite/CrsMatrix.h>
#include <PsimagLite/TypeToString.h>
#include <PsimagLite/Vector.h>

#include "Combinatorial.hh"

#include <bitset>
#include <climits>

namespace LanczosPlusPlus {

template <typename T1, typename T2>
std::ostream& operator<<(std::ostream& os, const std::pair<T1, T2>& p)
{
	os << p.first << " " << p.second;
	return os;
}

struct LanczosGlobals {

	typedef std::pair<int, int> PairIntType;
	typedef unsigned int long   WordType;

	static int FERMION_SIGN;

	enum
	{
		FERMION,
		BOSON
	};

	enum
	{
		SPIN_UP,
		SPIN_DOWN
	};

	enum ConnectionEnum
	{
		NONE
	};

	static int doSign(WordType a, SizeType i)
	{
		WordType mask = (1 << i) - 1;
		// Parity of single occupied between i and nsite-1
		return (PsimagLite::BitManip::count(a & mask) & 1) ? FERMION_SIGN : 1;
	}

	template <typename T> static void binRep(std::ostream& os, SizeType n, const T& a)
	{
		const char*        beg = reinterpret_cast<const char*>(&a);
		SizeType           len = sizeof(a);
		const char*        end = beg + len;
		PsimagLite::String buffer("");

		while (beg != end) {
			PsimagLite::String str(std::bitset<CHAR_BIT>(*(end - 1)).to_string());
			buffer += str;
			end--;
		}

		SizeType lmn = buffer.length();
		if (lmn >= n)
			lmn -= n;
		os << buffer.substr(lmn, n);
	}

	template <typename VectorWordType>
	static void printBasisBinary(std::ostream& os, SizeType n, VectorWordType& data)
	{
		for (SizeType i = 0; i < data.size(); i++) {
			binRep(os, n, data[i]);
			os << "\n";
		}

		os << "--------------\n";
	}

	static void printBasisBinary(std::ostream& os, SizeType n, SizeType h)
	{
		for (SizeType i = 0; i < h; ++i) {
			binRep(os, n, i);
			os << "\n";
		}

		os << "--------------\n";
	}

	template <typename WordType>
	static void
	printBasisDecimal(std::ostream& os, SizeType n, const std::vector<WordType>& data)
	{
		for (SizeType i = 0; i < data.size(); i++) {
			os << data[i] << " ";
			if (i > 0 && i % n == 0)
				std::cout << "\n";
		}

		os << "\n--------------\n";
	}

	template <typename WordType>
	static void printBasisDecimal(std::ostream&                                     os,
	                              SizeType                                          n,
	                              const std::vector<std::pair<WordType, WordType>>& data)
	{
		for (SizeType i = 0; i < data.size(); i++) {
			os << data[i].first << " " << data[i].second;
			if (i > 0 && i % n == 0)
				std::cout << "\n";
		}

		os << "\n--------------\n";
	}

	static void printBasisDecimal(std::ostream& os, SizeType n, SizeType h)
	{
		for (SizeType i = 0; i < h; ++i) {
			os << i << " ";
			if (i > 0 && i % n == 0)
				std::cout << "\n";
		}

		os << "\n--------------\n";
	}

	template <typename SomeVectorType>
	static typename PsimagLite::EnableIf<PsimagLite::IsVectorLike<SomeVectorType>::True,
	                                     void>::Type
	transform(SomeVectorType&                                                   gs,
	          SizeType                                                          offset,
	          SomeVectorType&                                                   gstmp,
	          const PsimagLite::CrsMatrix<typename SomeVectorType::value_type>& tr)
	{
		for (SizeType i = 0; i < gs.size(); i++) {
			assert(i + offset < gstmp.size());
			gstmp[i + offset] = gs[i];
		}

		PsimagLite::CrsMatrix<typename SomeVectorType::value_type> rT;
		transposeConjugate(rT, tr);
		gs.clear();
		gs.resize(tr.rows());
		multiply(gs, rT, gstmp);
	}

	static void doBitmask(SizeType total)
	{
		if (total == bitmask_.size())
			return;

		bitmask_.resize(total);
		bitmask_[0] = 1ul;
		for (SizeType i = 1; i < bitmask_.size(); i++)
			bitmask_[i] = bitmask_[i - 1] << 1;
	}

	static const WordType& bitmask(SizeType i)
	{
		assert(i < bitmask_.size());
		return bitmask_[i];
	}

	static void doCombinatorial(SizeType total)
	{
		if (total == comb_.size()) {
			return;
		}

		comb_.resize(total);
	}

	static const SizeType& combinatorial(SizeType i, SizeType j) { return comb_(i, j); }

private:

	static PsimagLite::Vector<WordType>::Type bitmask_;
	static Combinatorial                      comb_;
}; // LanczosGlobals

} // namespace LanczosPlusPlus
/*@}*/
#endif
