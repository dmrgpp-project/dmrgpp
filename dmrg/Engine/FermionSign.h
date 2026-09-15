// SPDX-FileCopyrightText: Copyright (c) 2009-2015, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file FermionSign.h
 *
 *  FIXME documentation
 */
#ifndef FEMION_SIGN_H
#define FEMION_SIGN_H
#include <PsimagLite/PackIndices.h>
#include <stdexcept>
#include <vector>

namespace Dmrg {

class FermionSign {

	using PackIndicesType = PsimagLite::PackIndices;
	using VectorSizeType  = PsimagLite::Vector<SizeType>::Type;
	using VectorBoolType  = PsimagLite::Vector<bool>::Type;

public:

	FermionSign(const VectorBoolType& signs)
	    : signs_(signs)
	{ }

	template <typename SomeBasisType>
	FermionSign(const SomeBasisType& basis, const VectorBoolType& signs)
	{
		if (basis.oldSigns().size() != basis.permutationInverse().size())
			err("FermionSign: Problem\n");

		SizeType n = basis.oldSigns().size();
		assert(basis.oldSigns().size() % signs.size() == 0);
		SizeType        nx = basis.oldSigns().size() / signs.size();
		PackIndicesType pack(nx);
		signs_.resize(nx);
		if (nx <= 1)
			return;
		for (SizeType x = 0; x < n; ++x) {
			SizeType x0 = 0;
			SizeType x1 = 0;
			pack.unpack(x0, x1, basis.permutation(x));
			assert(x1 < signs.size());
			bool parity1 = basis.oldSigns()[x];
			bool parity2 = signs[x1];
			assert(x0 < signs_.size());
			signs_[x0] = (parity1 != parity2);
		}
	}

	template <typename IoInputter>
	FermionSign(
	    IoInputter&        io,
	    PsimagLite::String prefix,
	    bool               bogus,
	    typename PsimagLite::EnableIf<PsimagLite::IsInputLike<IoInputter>::True, int>::Type = 0)
	{
		if (bogus)
			return;
		io.read(signs_, prefix);
	}

	int operator()(SizeType i, int f) const
	{
		assert(i < signs_.size());
		return (signs_[i]) ? f : 1;
	}

	template <typename IoOutputter>
	void
	write(IoOutputter&       io,
	      PsimagLite::String prefix,
	      typename PsimagLite::EnableIf<PsimagLite::IsOutputLike<IoOutputter>::True, int>::Type
	      = 0) const
	{
		io.write(signs_, prefix);
	}

	SizeType size() const { return signs_.size(); }

private:

	VectorBoolType signs_;
}; // class FermionSign
} // namespace Dmrg

/*@}*/
#endif // FEMION_SIGN_H
