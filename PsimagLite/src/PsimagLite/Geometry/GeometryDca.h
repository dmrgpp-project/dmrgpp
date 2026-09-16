// SPDX-FileCopyrightText: Copyright (c) 2014, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]

#ifndef LANCZOS_GEOMETRY_DCA_H
#define LANCZOS_GEOMETRY_DCA_H
#include <PsimagLite/Io/IoNg.h>
#include <PsimagLite/Vector.h>

namespace PsimagLite {

template <typename RealType_, typename GeometryType_> class GeometryDca {

public:

	using GeometryType = GeometryType_;
	using RealType     = RealType_;

	GeometryDca(const GeometryType& g, SizeType orbitals)
	    : enabled_(false)
	    , message_(g.term(0).label() + " with " + ttos(orbitals))
	{
		if (g.term(0).label() == "star" && orbitals == 4)
			enabled_ = true;
	}

	SizeType kSum(SizeType k1, SizeType k2) const
	{
		if (!enabled_)
			printErrorAndDie("kSum");
		return (k1 ^ k2);
	}

	SizeType kSustract(SizeType k1, SizeType k2) const
	{
		if (!enabled_)
			printErrorAndDie("kSustract");
		return (k1 ^ k2);
	}

	void write(String label1, IoNg::Out::Serializer& io) const
	{
		String label = label1 + "/GeometryDca";
		io.createGroup(label);
		io.write(label + "/enabled_", enabled_);
		io.write(label + "/message_", message_);
	}

private:

	void printErrorAndDie(String func) const
	{
		String str("GeometryDca:: Only valid for 2x2 clusters\n");
		str += "Not for " + message_ + " orbitals.\n";
		str += "\tWhile calling function " + func + "\n";
		throw RuntimeError(str);
	}

	bool   enabled_;
	String message_;
}; // class GeometryDca

} // namespace PsimagLite

#endif
