// SPDX-FileCopyrightText: Copyright (c) 2009-2014, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 2.0.0]
// [by G.A., Oak Ridge National Laboratory]


/** \ingroup DMRG */
/*@{*/

/*! \file ExtendedHubbard1Orb.h
 *  FIXME: Merge into Hubbard
 *  Hubbard + V_{ij} n_i n_j
 *
 */
#ifndef EXTENDED_HUBBARD_1ORB_H
#define EXTENDED_HUBBARD_1ORB_H
#include "../Models/HubbardOneBand/ModelHubbard.h"

namespace Dmrg {
//! Extended Hubbard for DMRG solver, uses ModelHubbard by containment
template <typename ModelBaseType> class ExtendedHubbard1Orb : public ModelBaseType {

public:

	using VectorSizeType          = typename ModelBaseType::VectorSizeType;
	using ModelHubbardType        = ModelHubbard<ModelBaseType>;
	using ModelHelperType         = typename ModelBaseType::ModelHelperType;
	using SuperGeometryType       = typename ModelBaseType::SuperGeometryType;
	using LeftRightSuperType      = typename ModelBaseType::LeftRightSuperType;
	using LinkType                = typename ModelBaseType::LinkType;
	using OperatorsType           = typename ModelHelperType::OperatorsType;
	using OperatorType            = typename OperatorsType::OperatorType;
	using RealType                = typename ModelHelperType::RealType;
	using QnType                  = typename ModelBaseType::QnType;
	using VectorQnType            = typename QnType::VectorQnType;
	using SparseMatrixType        = typename ModelHelperType::SparseMatrixType;
	using ComplexOrRealType       = typename SparseMatrixType::value_type;
	using MyBasis                 = typename ModelBaseType::MyBasis;
	using MyBasisWithOperators    = typename ModelBaseType::BasisWithOperatorsType;
	using HilbertBasisType        = typename ModelHubbardType::HilbertBasisType;
	using BlockType               = typename ModelHelperType::BlockType;
	using SolverParamsType        = typename ModelBaseType::SolverParamsType;
	using VectorType              = typename ModelBaseType::VectorType;
	using HilbertSpaceHubbardType = typename ModelHubbardType::HilbertSpaceHubbardType;
	using HilbertState            = typename HilbertSpaceHubbardType::HilbertState;
	using InputValidatorType      = typename ModelBaseType::InputValidatorType;
	using VectorOperatorType      = typename ModelBaseType::VectorOperatorType;
	using VectorHilbertStateType  = typename PsimagLite::Vector<HilbertState>::Type;
	using OpsLabelType            = typename ModelBaseType::OpsLabelType;
	using ModelTermType           = typename ModelBaseType::ModelTermType;
	using OpForLinkType           = typename ModelBaseType::OpForLinkType;

	ExtendedHubbard1Orb(const SolverParamsType&  solverParams,
	                    InputValidatorType&      io,
	                    SuperGeometryType const& geometry,
	                    PsimagLite::String       extension)
	    : ModelBaseType(solverParams, geometry, io)
	    , modelParameters_(io)
	    , superGeometry_(geometry)
	    , modelHubbard_(solverParams, io, geometry, extension)
	{ }

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const override
	{
		if (!io.doesGroupExist(label1))
			io.createGroup(label1);

		PsimagLite::String label = label1 + "/" + this->params().model;
		io.createGroup(label);
		modelParameters_.write(label, io);
		modelHubbard_.write(label, io);
	}

	void addDiagonalsInNaturalBasis(SparseMatrixType& hmatrix,
	                                const BlockType&  block,
	                                RealType          time) const override
	{
		modelHubbard_.addDiagonalsInNaturalBasis(hmatrix, block, time);
	}

	void fillLabeledOperators(VectorQnType& qns) override
	{
		modelHubbard_.fillLabeledOperators(qns);

		this->makeTrackable("n");
	}

	void fillModelLinks() override
	{
		modelHubbard_.fillModelLinks();

		ModelTermType& ninj = ModelBaseType::createTerm("ninj");
		OpForLinkType  n("n");

		ninj.push(n, 'N', n, 'N');
	}

private:

	//! Find n_i in the natural basis natBasis
	SparseMatrixType findOperatorMatrices(int i, const VectorHilbertStateType& natBasis) const
	{

		SizeType                                                  n = natBasis.size();
		PsimagLite::Matrix<typename SparseMatrixType::value_type> cm(n, n);

		for (SizeType ii = 0; ii < natBasis.size(); ii++) {
			HilbertState ket = natBasis[ii];
			cm(ii, ii)       = 0.0;
			for (SizeType sigma = 0; sigma < 2; sigma++)
				if (HilbertSpaceHubbardType::isNonZero(ket, i, sigma))
					cm(ii, ii) += 1.0;
		}

		SparseMatrixType creationMatrix(cm);
		return creationMatrix;
	}

	ParametersModelHubbard<RealType, QnType> modelParameters_;
	const SuperGeometryType&                 superGeometry_;
	ModelHubbardType                         modelHubbard_;
}; // class ExtendedHubbard1Orb

} // namespace Dmrg
/*@}*/
#endif // EXTENDED_HUBBARD_1ORB_H
