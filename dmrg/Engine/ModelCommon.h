// SPDX-FileCopyrightText: Copyright (c) 2009-2012-2018, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file ModelCommon.h
 *
 *  An abstract class to represent the strongly-correlated-electron models that
 *  can be used with the DmrgSolver
 *
 */

#ifndef MODEL_COMMON_H
#define MODEL_COMMON_H
#include <iostream>

#include "HamiltonianConnection.h"
#include "InputCheck.h"
#include "LabeledOperators.h"
#include "ModelLinks.h"
#include <PsimagLite/InputNg.h>
#include <PsimagLite/NoPthreads.h>
#include <PsimagLite/Profiling.h>
#include <PsimagLite/ProgressIndicator.h>
#include <PsimagLite/Sort.h>

namespace Dmrg {

template <typename ParametersType, typename SuperGeometryType, typename ModelHelperType>
class ModelCommon {

	using SparseMatrixType  = typename ModelHelperType::SparseMatrixType;
	using ComplexOrRealType = typename SparseMatrixType::value_type;
	using LinkType          = typename ModelHelperType::LinkType;

public:

	using VectorType             = typename PsimagLite::Vector<ComplexOrRealType>::Type;
	using InputValidatorType     = PsimagLite::InputNg<InputCheck>::Readable;
	using OperatorsType          = typename ModelHelperType::OperatorsType;
	using OperatorType           = typename OperatorsType::OperatorType;
	using Block                  = typename ModelHelperType::BlockType;
	using RealType               = typename ModelHelperType::RealType;
	using MyBasis                = typename ModelHelperType::BasisType;
	using BasisWithOperatorsType = typename ModelHelperType::BasisWithOperatorsType;
	using LabeledOperatorsType   = LabeledOperators<OperatorType>;
	using ModelLinksType         = ModelLinks<LabeledOperatorsType, SuperGeometryType>;
	using HamiltonianConnectionType
	    = HamiltonianConnection<ModelLinksType, ModelHelperType, ParametersType>;
	using VectorLinkType       = typename HamiltonianConnectionType::VectorLinkType;
	using LeftRightSuperType   = typename ModelHelperType::LeftRightSuperType;
	using VectorOperatorType   = typename PsimagLite::Vector<OperatorType>::Type;
	using VectorVectorLinkType = typename PsimagLite::Vector<VectorLinkType>::Type;
	using VectorSizeType       = typename HamiltonianConnectionType::VectorSizeType;
	using VerySparseMatrixType = typename HamiltonianConnectionType::VerySparseMatrixType;

	ModelCommon(const ParametersType& params, const SuperGeometryType& superGeometry)
	    : params_(params)
	    , superGeometry_(superGeometry)
	    , progress_("ModelCommon")
	{
		if (params.options.isSet("OperatorsChangeAll"))
			OperatorsType::setChangeAll(true);
	}

	const ParametersType& params() const { return params_; }

	const SuperGeometryType& superGeometry() const { return superGeometry_; }

	void addConnectionsInNaturalBasis(SparseMatrixType&         hmatrix,
	                                  const VectorOperatorType& cm,
	                                  const Block&              block,
	                                  RealType                  time) const
	{
		if (block.size() != 1)
			err("addConnectionsInNaturalBasis(): unimplemented\n");
	}

	const ParametersType&         params_;
	const SuperGeometryType&      superGeometry_;
	PsimagLite::ProgressIndicator progress_;
}; // class ModelCommon
} // namespace Dmrg
/*@}*/
#endif
