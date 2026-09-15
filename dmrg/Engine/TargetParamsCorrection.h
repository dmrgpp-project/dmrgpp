// SPDX-FileCopyrightText: Copyright (c) 2009, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 2.0.0]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file TargetParamsCorrection.h
 *
 * Contains the parameter a of PRB 72 180403(R) (2005)
 *
 */
#ifndef TARGET_PARAMS_CORRECTION_H
#define TARGET_PARAMS_CORRECTION_H

#include "TargetParamsBase.h"
#include <vector>

namespace Dmrg {
// Coordinates reading of TargetSTructure from input file
template <typename ModelType> class TargetParamsCorrection : public TargetParamsBase<ModelType> {

	using BaseType = TargetParamsBase<ModelType>;

public:

	using RealType = typename BaseType::RealType;

	template <typename IoInputter>
	TargetParamsCorrection(IoInputter& io, PsimagLite::String targeting)
	    : BaseType(targeting)
	{
		io.readline(correctionA_, "CorrectionA=");
	}

	virtual SizeType memResolv(PsimagLite::MemResolv&, SizeType, PsimagLite::String = "") const
	{
		return 0;
	}

	SizeType sites() const override { return 0; }

	RealType correctionA() const override { return correctionA_; }

	SizeType sectorIndex() const override
	{
		throw PsimagLite::RuntimeError("sectorIndex called for gs\n");
	}

	SizeType levelIndex() const override
	{
		throw PsimagLite::RuntimeError("levelIndex called for gs\n");
	}

	void write(PsimagLite::String label, PsimagLite::IoSerializer& ioSerializer) const
	{
		ioSerializer.createGroup(label);
		ioSerializer.write(label + "/correctionA_", correctionA_);
	}

private:

	RealType correctionA_;
}; // class TargetParamsCorrection

template <typename ModelType>
inline std::ostream& operator<<(std::ostream& os, const TargetParamsCorrection<ModelType>& t)
{
	os << "TargetParams.type=correction\n";
	os << "TargetCorrection.correctionA=" << t.correctionA() << "\n";
	return os;
}
} // namespace Dmrg

/*@}*/
#endif // TARGET_PARAMS_CORRECTION_H
