// SPDX-FileCopyrightText: Copyright (c) 2009-2013, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file SuperHubbardExtended.h
 *
 *  Hubbard + V_{ij} n_i n_j
 *
 */
#ifndef SUPER_HUBBARD_EXTENDED_H
#define SUPER_HUBBARD_EXTENDED_H
#include "../Models/ExtendedHubbard1Orb/ExtendedHubbard1Orb.h"

namespace Dmrg {
//! Extended Hubbard for DMRG solver, uses ExtendedHubbard1Orb by containment
template <typename ModelBaseType> class SuperHubbardExtended : public ModelBaseType {

public:

	using VectorSizeType          = typename ModelBaseType::VectorSizeType;
	using ModelHubbardType        = ExtendedHubbard1Orb<ModelBaseType>;
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
	using BasisType               = typename ModelBaseType::MyBasis;
	using MyBasisWithOperators    = typename ModelBaseType::BasisWithOperatorsType;
	using HilbertBasisType        = typename ModelHubbardType::HilbertBasisType;
	using BlockType               = typename ModelHelperType::BlockType;
	using SolverParamsType        = typename ModelBaseType::SolverParamsType;
	using VectorType              = typename ModelBaseType::VectorType;
	using HilbertSpaceHubbardType = typename ModelHubbardType::HilbertSpaceHubbardType;
	using HilbertState            = typename HilbertSpaceHubbardType::HilbertState;
	using InputValidatorType      = typename ModelBaseType::InputValidatorType;
	using VectorOperatorType      = typename ModelBaseType::VectorOperatorType;
	using ModelTermType           = typename ModelBaseType::ModelTermType;
	using VectorHilbertStateType  = typename PsimagLite::Vector<HilbertState>::Type;
	using VectorSizeTypeType      = PsimagLite::Vector<SizeType>::Type;
	using OpsLabelType            = typename ModelBaseType::OpsLabelType;
	using OpForLinkType           = typename ModelBaseType::OpForLinkType;

	SuperHubbardExtended(const SolverParamsType&  solverParams,
	                     InputValidatorType&      io,
	                     const SuperGeometryType& geometry)
	    : ModelBaseType(solverParams, geometry, io)
	    , modelParameters_(io)
	    , extendedHubbard_(solverParams, io, geometry, "")
	{ }

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const override
	{
		if (!io.doesGroupExist(label1))
			io.createGroup(label1);

		PsimagLite::String label = label1 + "/" + this->params().model;
		io.createGroup(label);
		modelParameters_.write(label, io);
		extendedHubbard_.write(label, io);
	}

	void addDiagonalsInNaturalBasis(SparseMatrixType& hmatrix,
	                                const BlockType&  block,
	                                RealType          time) const override
	{
		extendedHubbard_.addDiagonalsInNaturalBasis(hmatrix, block, time);
	}

protected:

	void fillLabeledOperators(VectorQnType& qns) override
	{
		extendedHubbard_.fillLabeledOperators(qns);

		this->makeTrackable("splus");
		this->makeTrackable("sz");
	}

	void fillModelLinks() override
	{
		extendedHubbard_.fillModelLinks();

		const bool isSu2 = BasisType::useSu2Symmetry();

		ModelTermType& spsm = ModelBaseType::createTerm("SplusSminus");

		OpForLinkType splus("splus");

		auto valueModiferTerm0
		    = [isSu2](ComplexOrRealType& value) { value *= (isSu2) ? -0.5 : 0.5; };

		spsm.push(splus, 'N', splus, 'C', valueModiferTerm0);

		ModelTermType& szsz = ModelBaseType::createTerm("szsz");

		if (!isSu2) {
			OpForLinkType sz("sz");
			szsz.push(sz, 'N', sz, 'N', typename ModelTermType::Su2Properties(2, 0.5));
		} else {
			auto valueModifierTermOther = [isSu2](ComplexOrRealType& value)
			{
				if (isSu2)
					value = -value;
			};
			spsm.push(splus, 'N', splus, 'C', valueModifierTermOther);
		}
	}

private:

	//! Full hamiltonian from creation matrices cm
	void addSiSj(SparseMatrixType&, const VectorOperatorType&, const BlockType& block) const
	{
		// Assume block.size()==1 and then problem solved!!
		// there are no connection if there's only one site ;-)
		assert(block.size() == 1);
	}

	ParametersModelHubbard<RealType, QnType> modelParameters_;
	ModelHubbardType                         extendedHubbard_;
}; // class SuperHubbardExtended

} // namespace Dmrg
/*@}*/
#endif // SUPER_HUBBARD_EXTENDED_H
