// SPDX-FileCopyrightText: Copyright (c) 2009, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 2.0.0]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file TargetParamsDynamic.h
 *
 *  This is a structure to represent the parameters of the
 *  Dynamic DMRG algorithm.
 *  Don't add functions to this class because
 *  this class's data is all public
 */
#ifndef TARGET_PARAMS_DYNAMIC_H
#define TARGET_PARAMS_DYNAMIC_H

#include "TargetParamsCommon.h"

namespace Dmrg {
// Coordinates reading of TargetSTructure from input file
template <typename ModelType> class TargetParamsDynamic : public TargetParamsCommon<ModelType> {
public:

	using TargetParamsCommonType = TargetParamsCommon<ModelType>;
	using RealType               = typename ModelType::RealType;

	using OperatorType     = typename ModelType::OperatorType;
	using PairType         = typename OperatorType::PairType;
	using SparseMatrixType = typename OperatorType::StorageType;
	using ComplexOrReal    = typename SparseMatrixType::value_type;
	using MatrixType       = PsimagLite::Matrix<ComplexOrReal>;

	template <typename IoInputter>
	TargetParamsDynamic(IoInputter& io, PsimagLite::String targeting, const ModelType& model)
	    : TargetParamsCommonType(io, targeting, model)
	{
		io.readline(type_, "DynamicDmrgType=");
	}

	virtual SizeType memResolv(PsimagLite::MemResolv&, SizeType, PsimagLite::String = "") const
	{
		return 0;
	}

	SizeType type() const override { return type_; }

private:

	SizeType type_;

}; // class TargetParamsDynamic

template <typename ModelType>
inline std::ostream& operator<<(std::ostream& os, const TargetParamsDynamic<ModelType>& t)
{
	os << "TargetParams.type=DynamicDmrg\n";
	const typename TargetParamsDynamic<ModelType>::TargetParamsCommonType& tp = t;
	os << tp;
	os << "DynamicDmrgType=" << t.type() << "\n";
	return os;
}
} // namespace Dmrg

/*@}*/
#endif // TARGET_PARAMS_DYNAMIC_H
