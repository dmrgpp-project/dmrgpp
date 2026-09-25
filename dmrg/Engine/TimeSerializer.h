// SPDX-FileCopyrightText: Copyright (c) 2009, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 2.0.0]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file TimeSerializer.h
 *
 *  Serialize time data
 */
#ifndef TIME_SERIAL_H
#define TIME_SERIAL_H

#include "StageEnum.h"
#include <PsimagLite/Io/IoSelector.h>
#include <PsimagLite/TypeToString.h>
#include <PsimagLite/Vector.h>

namespace Dmrg {

template <typename VectorType> class TimeSerializer {

public:

	using VectorElementType   = typename VectorType::value_type;
	using RealType            = typename PsimagLite::Real<VectorElementType>::Type;
	using VectorStageEnumType = typename PsimagLite::Vector<StageEnum>::Type;
	using VectorVectorType    = typename PsimagLite::Vector<VectorType*>::Type;

	template <typename SomeAoeType>
	TimeSerializer(SizeType           currentTimeStep,
	               RealType           currentTime,
	               SizeType           site,
	               const SomeAoeType& aoe,
	               PsimagLite::String name)
	    : currentTimeStep_(currentTimeStep)
	    , currentTime_(currentTime)
	    , site_(site)
	    , targetVectors_(aoe.tvs())
	    , stages_(aoe.stages())
	    , name_(name)
	    , owner_(false)
	{
		const SizeType n = targetVectors_.size();
		for (SizeType i = 0; i < n; ++i)
			targetVectors_[i] = const_cast<VectorType*>(&aoe.targetVectors(i));
	}

	TimeSerializer(typename PsimagLite::IoSelector::In& io, PsimagLite::String prefix)
	    : owner_(true)
	{
		prefix += "/TimeSerializer/";

		PsimagLite::String s = prefix + "CurrentTimeStep";
		io.read(currentTimeStep_, s);

		s = prefix + "Time";
		io.read(currentTime_, s);

		s      = prefix + "TargetCentralSite";
		int xi = 0;
		io.read(xi, s);
		if (xi < 0)
			err("TimeSerializer:: site cannot be negative\n");
		site_ = xi;

		s = prefix + "TNUMBEROFVECTORS";
		io.read(xi, s);
		if (xi <= 0)
			err("TimeSerializer:: n. of vectors must be positive\n");
		targetVectors_.clear();
		for (int i = 0; i < xi; ++i) {
			VectorType* v = new VectorType();
			s             = prefix + "targetVector" + ttos(i);
			v->read(io, s);
			targetVectors_.push_back(v);
		}

		s = prefix + "Stages";
		io.read(stages_, s);

		try {
			io.read(name_, prefix + "Name");
		} catch (...) {
			// reading an old file, set name to LEGACY
			name_ = "LEGACY";
		}
	}

	~TimeSerializer()
	{
		if (!owner_)
			return;
		const SizeType n = targetVectors_.size();
		for (SizeType i = 0; i < n; ++i) {
			delete targetVectors_[i];
			targetVectors_[i] = nullptr;
		}
	}

	void write(PsimagLite::IoSelector::Out& io, PsimagLite::String prefix) const
	{
		prefix += "/TimeSerializer";
		io.createGroup(prefix);
		prefix += "/";

		io.write(currentTime_, prefix + "Time");
		io.write(currentTimeStep_, prefix + "CurrentTimeStep");
		io.write(site_, prefix + "TargetCentralSite");
		io.write(targetVectors_.size(), prefix + "TNUMBEROFVECTORS");

		for (SizeType i = 0; i < targetVectors_.size(); i++) {
			PsimagLite::String label = "targetVector" + ttos(i);
			targetVectors_[i]->write(io, prefix + label);
		}

		io.write(stages_, prefix + "Stages");
		io.write(name_, prefix + "Name");
	}

	SizeType numberOfVectors() const { return targetVectors_.size(); }

	SizeType currentTimeStep() const { return currentTimeStep_; }

	RealType time() const { return currentTime_; }

	SizeType site() const { return site_; }

	PsimagLite::String name() const { return name_; }

	const VectorType& vector(SizeType i) const
	{
		if (i >= targetVectors_.size())
			err("TimeSerializer: Not so many time vectors\n");

		if (!targetVectors_[i])
			err("TimeSerializer: FATAL\n");

		return *targetVectors_[i];
	}

	const VectorStageEnumType& stages() const { return stages_; }

private:

	SizeType            currentTimeStep_;
	RealType            currentTime_;
	SizeType            site_;
	VectorVectorType    targetVectors_;
	VectorStageEnumType stages_;
	PsimagLite::String  name_;
	bool                owner_;
}; // class TimeSerializer
} // namespace Dmrg

/*@}*/
#endif
