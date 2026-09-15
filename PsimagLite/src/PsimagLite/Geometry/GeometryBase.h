// SPDX-FileCopyrightText: Copyright (c) 2009-2013, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup PsimagLite */
/*@{*/

/*! \file GeometryBase.h
 *
 *  Well, I need to read that chapter in
 *  Alexandrescu's "Modern C++ design" again to have
 *  a decent factory here, but this will have to do for now
 *
 */
#ifndef GEOMETRY_BASE_H
#define GEOMETRY_BASE_H

#include <PsimagLite/InputNg.h>

namespace PsimagLite {

template <typename ComplexOrRealType, typename InputType> class GeometryBase {

	using PairType   = std::pair<SizeType, SizeType>;
	using MatrixType = Matrix<ComplexOrRealType>;

public:

	virtual ~GeometryBase() { }

	template <class Archive> void write(Archive&, const unsigned int) { }

	virtual SizeType dirs() const = 0;

	virtual SizeType handle(SizeType i, SizeType j) const = 0;

	virtual SizeType getVectorSize(SizeType dirId) const = 0;

	virtual bool connected(SizeType i1, SizeType i2) const = 0;

	virtual SizeType calcDir(SizeType i1, SizeType i2) const = 0;

	virtual bool fringe(SizeType i, SizeType smax, SizeType emin) const = 0;

	virtual SizeType getSubstituteSite(SizeType smax, SizeType emin, SizeType siteNew2) const
	    = 0;

	virtual String label() const = 0;

	virtual SizeType length(SizeType i) const = 0;

	virtual SizeType translate(SizeType site, SizeType dir, SizeType amount) const = 0;

	virtual SizeType maxConnections() const = 0;

	virtual SizeType findReflection(SizeType site) const = 0;

	virtual void set(MatrixType&, SizeType) const
	{
		throw RuntimeError("GeometryBase::set() unimplemented for derived class\n");
	}

	virtual int index(SizeType i1, SizeType edof1, SizeType edofTotal) const
	{
		assert(edof1 < edofTotal);
		return edof1 + i1 * edofTotal;
	}

	virtual SizeType matrixRank(SizeType linSize, SizeType maxEdof) const
	{
		return linSize * maxEdof;
	}

	virtual int signChange(SizeType, SizeType) const { return 1; }

	virtual SizeType orbitals(SizeType orbs, SizeType) const { return orbs; }

protected:

	SizeType unimplemented(const String& str) const
	{
		String str2 = "unimplemented " + str + "\n";
		throw RuntimeError(str2);
	}

	bool neighbors(SizeType i1, SizeType i2, bool periodic, SizeType period) const
	{
		SizeType imin = (i1 < i2) ? i1 : i2;
		SizeType imax = (i1 > i2) ? i1 : i2;
		bool     b    = (imax - imin == 1);
		if (!periodic)
			return b;
		bool b2 = (imax - imin == period);
		return (b || b2);
	}
}; // class GeometryBase
} // namespace PsimagLite

/*@}*/
#endif // GEOMETRY_BASE_H
