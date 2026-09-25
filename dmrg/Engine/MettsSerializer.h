// SPDX-FileCopyrightText: Copyright (c) 2009-2014, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file MettsSerializer.h
 *
 *  Serialize time data
 */
#ifndef METTS_SERIALIZER_H
#define METTS_SERIALIZER_H

#include <PsimagLite/Io/IoSelector.h>
#include <PsimagLite/TypeToString.h>

namespace Dmrg {

template <typename VectorType> class MettsSerializer {

	using VectorElementType = typename VectorType::value_type;
	using RealType          = typename PsimagLite::Real<VectorElementType>::Type;

public:

	// Unfortunately we need a default ctor
	// to build an array of these
	MettsSerializer() { }

	MettsSerializer(RealType                                             currentBeta,
	                SizeType                                             site,
	                const typename PsimagLite::Vector<VectorType>::Type& targetVectors)
	    : currentBeta_(currentBeta)
	    , site_(site)
	    , targetVectors_(targetVectors)
	{ }

	MettsSerializer(typename PsimagLite::IoSelector::In& io)
	{
		RealType           x = 0;
		PsimagLite::String s = "BETA";

		io.read(x, s);

		if (x < 0)
			throw PsimagLite::RuntimeError(
			    "MettsSerializer:: time cannot be negative\n");

		currentBeta_ = x;

		s      = "TargetCentralSite";
		int xi = 0;
		io.read(xi, s);
		if (xi < 0)
			throw PsimagLite::RuntimeError(
			    "MettsSerializer:: site cannot be negative\n");

		site_ = xi;

		s = "TNUMBEROFVECTORS";
		io.read(xi, s);
		if (xi <= 0)
			throw PsimagLite::RuntimeError(
			    "MettsSerializer:: n. of vectors must be positive\n");

		targetVectors_.resize(xi);
		for (SizeType i = 0; i < targetVectors_.size(); i++) {
			s = "targetVector" + ttos(i);
			targetVectors_[i].read(io, s);
		}
	}

	SizeType size(SizeType i = 0) const { return targetVectors_[i].size(); }

	RealType beta() const { return currentBeta_; }

	SizeType site() const { return site_; }

	const VectorType& vector(SizeType i = 0) const { return targetVectors_[i]; }

	template <typename IoOutputter>
	void
	write(IoOutputter& io,
	      typename PsimagLite::EnableIf<PsimagLite::IsOutputLike<IoOutputter>::True, int>::Type
	      = 0) const
	{
		PsimagLite::String s = "BETA=" + ttos(currentBeta_);
		io.printline(s);
		s = "TargetCentralSite=" + ttos(site_);
		io.printline(s);
		s = "TNUMBEROFVECTORS=" + ttos(targetVectors_.size());
		io.printline(s);
		for (SizeType i = 0; i < targetVectors_.size(); i++) {
			PsimagLite::String label
			    = "targetVector" + ttos(i) + "_" + ttos(currentBeta_);
			targetVectors_[i].write(io, label);
		}
	}

private:

	RealType                                      currentBeta_;
	SizeType                                      site_;
	typename PsimagLite::Vector<VectorType>::Type targetVectors_;
}; // class MettsSerializer
} // namespace Dmrg

/*@}*/
#endif // METTS_SERIALIZER_H
