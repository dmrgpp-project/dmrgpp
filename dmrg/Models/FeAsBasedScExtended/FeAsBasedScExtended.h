// SPDX-FileCopyrightText: Copyright (c) 2009, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 2.0.0]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file FeBasedScExtedned.h
 *
 *  An implementation of a Hubbard model for Fe-based superconductors
 *  to use with the DmrgSolver
 *  This extends the FeAsBasedSc model to include JNN and JNNN couplings
 *  FIXME: Merge into FeAsBasedSc
 *
 */
#ifndef FEAS_BASED_SC_EX
#define FEAS_BASED_SC_EX
#include "../Models/FeAsModel/ModelFeBasedSc.h"

namespace Dmrg {

template <typename ModelBaseType> class FeAsBasedScExtended : public ModelBaseType {

public:

	using VectorSizeType       = typename ModelBaseType::VectorSizeType;
	using ModelFeAsType        = ModelFeBasedSc<ModelBaseType>;
	using HilbertState         = typename ModelFeAsType::HilbertState;
	using HilbertBasisType     = typename ModelFeAsType::HilbertBasisType;
	using ModelHelperType      = typename ModelBaseType::ModelHelperType;
	using SuperGeometryType    = typename ModelBaseType::SuperGeometryType;
	using LeftRightSuperType   = typename ModelBaseType::LeftRightSuperType;
	using LinkType             = typename ModelBaseType::LinkType;
	using OperatorsType        = typename ModelHelperType::OperatorsType;
	using OperatorType         = typename OperatorsType::OperatorType;
	using VectorOperatorType   = typename PsimagLite::Vector<OperatorType>::Type;
	using RealType             = typename ModelHelperType::RealType;
	using QnType               = typename ModelBaseType::QnType;
	using VectorQnType         = typename QnType::VectorQnType;
	using SparseMatrixType     = typename ModelHelperType::SparseMatrixType;
	using ComplexOrRealType    = typename SparseMatrixType::value_type;
	using Su2RelatedType       = typename OperatorType::Su2RelatedType;
	using BasisType            = typename ModelBaseType::MyBasis;
	using MyBasisWithOperators = typename ModelBaseType::BasisWithOperatorsType;
	using BlockType            = typename BasisType::BlockType;
	using SolverParamsType     = typename ModelBaseType::SolverParamsType;
	using VectorType           = typename ModelBaseType::VectorType;
	using InputValidatorType   = typename ModelBaseType::InputValidatorType;
	using OpsLabelType         = typename ModelBaseType::OpsLabelType;
	using OpForLinkType        = typename ModelBaseType::OpForLinkType;
	using ModelTermType        = typename ModelBaseType::ModelTermType;

	static const SizeType SPIN_UP   = ModelFeAsType::SPIN_UP;
	static const SizeType SPIN_DOWN = ModelFeAsType::SPIN_DOWN;

	FeAsBasedScExtended(const SolverParamsType&  solverParams,
	                    InputValidatorType&      io,
	                    const SuperGeometryType& geometry)
	    : ModelBaseType(solverParams, geometry, io)
	    , modelParameters_(io)
	    , modelFeAs_(solverParams, io, geometry, "")
	    , orbitals_(modelParameters_.orbitals)
	{ }

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const override
	{
		if (!io.doesGroupExist(label1))
			io.createGroup(label1);

		PsimagLite::String label = label1 + "/" + this->params().model;
		io.createGroup(label);
		modelParameters_.write(label, io);
		modelFeAs_.write(label, io);
		io.write(label + "/orbitals_", orbitals_);
	}

	void addDiagonalsInNaturalBasis(SparseMatrixType& hmatrix,
	                                const BlockType&  block,
	                                RealType          time) const override
	{
		modelFeAs_.addDiagonalsInNaturalBasis(hmatrix, block, time);
	}

protected:

	void fillLabeledOperators(VectorQnType& qns) override
	{
		modelFeAs_.fillLabeledOperators(qns);
		SizeType                                        site = 0;
		BlockType                                       block(1, site);
		typename PsimagLite::Vector<OperatorType>::Type creationMatrix;
		setOperatorMatricesInternal(creationMatrix, qns, block);

		OpsLabelType& splus = this->createOpsLabel("naturalSplus");
		SizeType      x     = 2 * orbitals_;
		assert(x < creationMatrix.size());
		splus.push(creationMatrix[x]);

		OpsLabelType& sminus = this->createOpsLabel("naturalSminus");
		x                    = 2 * orbitals_;
		assert(x < creationMatrix.size());
		creationMatrix[x].dagger();
		sminus.push(creationMatrix[x]);
		creationMatrix[x].dagger();

		OpsLabelType& sz = this->createOpsLabel("naturalSz");
		x                = 2 * orbitals_ + 1;
		assert(x < creationMatrix.size());
		sz.push(creationMatrix[x]);

		this->makeTrackable("naturalSplus");
		this->makeTrackable("naturalSz");
	}

	void fillModelLinks() override
	{
		modelFeAs_.fillModelLinks();

		bool isSu2 = BasisType::useSu2Symmetry();

		ModelTermType& spsm = ModelBaseType::createTerm("SplusSminus");

		OpForLinkType splus("naturalSplus");

		auto valueModiferTerm0
		    = [isSu2](ComplexOrRealType& value) { value *= (isSu2) ? -0.5 : 0.5; };

		spsm.push(splus, 'N', splus, 'C', valueModiferTerm0);

		ModelTermType& szsz = ModelBaseType::createTerm("szsz");

		OpForLinkType sz("naturalSz");
		szsz.push(sz, 'N', sz, 'N');
	}

private:

	//! set creation matrices for sites in block
	void setOperatorMatricesInternal(VectorOperatorType& creationMatrix,
	                                 VectorQnType&       qns,
	                                 const BlockType&    block) const
	{
		blockIsSize1OrThrow(block);

		modelFeAs_.setQns(qns);
		modelFeAs_.setOperatorMatricesInternal(creationMatrix, block);

		// add S^+_i to creationMatrix
		setSplus(creationMatrix, block);

		// add S^z_i to creationMatrix
		setSz(creationMatrix, block);
	}

	// add S^+_i to creationMatrix
	void setSplus(typename PsimagLite::Vector<OperatorType>::Type& creationMatrix,
	              const BlockType&                                 block) const
	{
		SparseMatrixType m;
		cDaggerC(m, creationMatrix, block, 1.0, SPIN_UP, SPIN_DOWN);
		Su2RelatedType su2related;
		SizeType       offset = 2 * orbitals_;
		su2related.source.push_back(offset);
		su2related.source.push_back(offset + 1);
		su2related.source.push_back(offset);
		su2related.transpose.push_back(-1);
		su2related.transpose.push_back(-1);
		su2related.transpose.push_back(1);
		su2related.offset = 1;

		OperatorType sPlus(m,
		                   ProgramGlobals::FermionOrBosonEnum::BOSON,
		                   typename OperatorType::PairType(2, 2),
		                   -1,
		                   su2related);
		creationMatrix.push_back(sPlus);
	}

	// add S^z_i to creationMatrix
	void setSz(typename PsimagLite::Vector<OperatorType>::Type& creationMatrix,
	           const BlockType&                                 block) const
	{
		SparseMatrixType m1, m2;
		cDaggerC(m1, creationMatrix, block, 0.5, SPIN_UP, SPIN_UP);
		cDaggerC(m2, creationMatrix, block, -0.5, SPIN_DOWN, SPIN_DOWN);
		Su2RelatedType   su2related2;
		SparseMatrixType m = m1;
		m += m2;
		OperatorType sz(m,
		                ProgramGlobals::FermionOrBosonEnum::BOSON,
		                typename OperatorType::PairType(2, 1),
		                1.0 / sqrt(2.0),
		                su2related2);
		creationMatrix.push_back(sz);
	}

	// add S^+_i to creationMatrix
	void cDaggerC(SparseMatrixType&                                      sum,
	              const typename PsimagLite::Vector<OperatorType>::Type& creationMatrix,
	              const BlockType&,
	              RealType value,
	              SizeType spin1,
	              SizeType spin2) const
	{
		SparseMatrixType tmpMatrix, tmpMatrix2;
		for (SizeType orbital = 0; orbital < orbitals_; orbital++) {
			transposeConjugate(tmpMatrix2,
			                   creationMatrix[orbital + spin2 * orbitals_].getCRS());
			multiply(tmpMatrix,
			         creationMatrix[orbital + spin1 * orbitals_].getCRS(),
			         tmpMatrix2);

			if (orbital == 0)
				sum = value * tmpMatrix;
			else
				sum += value * tmpMatrix;
		}
	}

	// add J_{ij} S^+_i S^-_j + S^-_i S^+_j to Hamiltonia
	void addSplusSminus(SparseMatrixType&,
	                    const typename PsimagLite::Vector<OperatorType>::Type&,
	                    const BlockType&) const
	{
		// nothing if block.size == 1
	}

	// add J_{ij} S^z_i S^z_j to Hamiltonian
	void addSzSz(SparseMatrixType&,
	             const typename PsimagLite::Vector<OperatorType>::Type&,
	             const BlockType&) const
	{
		// nothing if block.size == 1
	}

	void blockIsSize1OrThrow(const BlockType& block) const
	{
		if (block.size() == 1)
			return;
		throw PsimagLite::RuntimeError("FeAsBasedExtended:: blocks must be of size 1\n");
	}

	ParametersModelFeAs<RealType, QnType> modelParameters_;
	ModelFeAsType                         modelFeAs_;
	SizeType                              orbitals_;
}; // class FeAsBasedScExtended

} // namespace Dmrg
/*@}*/
#endif // FEAS_BASED_SC_EX
