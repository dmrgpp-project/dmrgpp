// SPDX-FileCopyrightText: Copyright (c) 2009, 2017-2019, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file GaugeSpin.h
 *
 *  An implementation of the Quantum Heisenberg Model to use with  DmrgSolver
 *
 */

#ifndef DMRG_MODEL_GAUGESPIN_H
#define DMRG_MODEL_GAUGESPIN_H

#include "ParametersGaugeSpin.h"
#include "ProgramGlobals.h"
#include "SuperOpHelperPlaquette.h"
#include "Utils.h"
#include "VerySparseMatrix.h"
#include <PsimagLite/CrsMatrix.h>
#include <algorithm>

namespace Dmrg {

template <typename ModelBaseType> class GaugeSpin : public ModelBaseType {

public:

	using ModelHelperType       = typename ModelBaseType::ModelHelperType;
	using BasisType             = typename ModelHelperType::BasisType;
	using SuperGeometryType     = typename ModelBaseType::SuperGeometryType;
	using LeftRightSuperType    = typename ModelBaseType::LeftRightSuperType;
	using LinkType              = typename ModelBaseType::LinkType;
	using OperatorsType         = typename ModelHelperType::OperatorsType;
	using RealType              = typename ModelHelperType::RealType;
	using VectorType            = typename ModelBaseType::VectorType;
	using QnType                = typename ModelBaseType::QnType;
	using VectorQnType          = typename ModelBaseType::VectorQnType;
	using BlockType             = typename ModelBaseType::BlockType;
	using SolverParamsType      = typename ModelBaseType::SolverParamsType;
	using SparseMatrixType      = typename ModelHelperType::SparseMatrixType;
	using ComplexOrRealType     = typename SparseMatrixType::value_type;
	using WordType              = unsigned int long;
	using InputValidatorType    = typename ModelBaseType::InputValidatorType;
	using MatrixType            = PsimagLite::Matrix<ComplexOrRealType>;
	using VectorSizeType        = typename PsimagLite::Vector<SizeType>::Type;
	using VectorRealType        = typename ModelBaseType::VectorRealType;
	using ModelTermType         = typename ModelBaseType::ModelTermType;
	using HilbertBasisType      = typename PsimagLite::Vector<SizeType>::Type;
	using OperatorType          = typename OperatorsType::OperatorType;
	using PairType              = typename OperatorType::PairType;
	using VectorOperatorType    = typename PsimagLite::Vector<OperatorType>::Type;
	using MyBasis               = typename ModelBaseType::MyBasis;
	using MyBasisWithOperators  = typename ModelBaseType::BasisWithOperatorsType;
	using OpsLabelType          = typename ModelBaseType::OpsLabelType;
	using OpForLinkType         = typename ModelBaseType::OpForLinkType;
	using SuperOpHelperBaseType = typename ModelBaseType::SuperOpHelperBaseType;
	using SuperOpHelperPlaquetteType
	    = SuperOpHelperPlaquette<SuperGeometryType, SolverParamsType>;

	static const SizeType TWICE_THE_SPIN = 2;

	GaugeSpin(const SolverParamsType&  solverParams,
	          InputValidatorType&      io,
	          const SuperGeometryType& geometry)
	    : ModelBaseType(solverParams, geometry, io)
	    , modelParameters_(io)
	    , superOpHelperPlaquette_(nullptr)
	{ }

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const override
	{
		if (!io.doesGroupExist(label1))
			io.createGroup(label1);

		PsimagLite::String label = label1 + "/" + this->params().model;
		io.createGroup(label);
		modelParameters_.write(label, io);
	}

	void addDiagonalsInNaturalBasis(SparseMatrixType& hmatrix,
	                                const BlockType&  block,
	                                RealType          time) const override
	{
		ModelBaseType::additionalOnSiteHamiltonian(hmatrix, block, time);

		SizeType linSize = ModelBaseType::superGeometry().numberOfSites();
		assert(block.size() == 1);

		SizeType site = block[0];

		const OperatorType& sz = ModelBaseType::naturalOperator("sz", site, 0);

		if (modelParameters_.magneticFieldV.size() == linSize) {

			RealType tmp = modelParameters_.magneticFieldV[site];
			hmatrix += tmp * sz.getCRS();
		}
	}

protected:

	void fillLabeledOperators(VectorQnType& qns) override
	{
		const SizeType   total = TWICE_THE_SPIN + 1;
		HilbertBasisType natBasis(total);
		for (SizeType i = 0; i < total; ++i)
			natBasis[i] = i;

		setSymmetryRelated(qns, natBasis);

		// Set the operators S^+_i in the natural basis
		SparseMatrixType tmpMatrix = findSplusMatrices(0, natBasis);

		typename OperatorType::Su2RelatedType su2related;

		OperatorType myOp(tmpMatrix,
		                  ProgramGlobals::FermionOrBosonEnum::BOSON,
		                  PairType(2, 2),
		                  -1,
		                  su2related);
		this->createOpsLabel("splus").push(myOp);
		// this->makeTrackable("splus");

		myOp.dagger();
		this->createOpsLabel("sminus").push(myOp);

		// Set the operators S^z_i in the natural basis
		tmpMatrix = findSzMatrices(0, natBasis);
		typename OperatorType::Su2RelatedType su2related2;
		OperatorType                          myOp2(tmpMatrix,
                                   ProgramGlobals::FermionOrBosonEnum::BOSON,
                                   PairType(2, 1),
                                   1.0 / sqrt(2.0),
                                   su2related2);
		this->createOpsLabel("sz").push(myOp2);
		// this->makeTrackable("sz");

		// Set the operators S^x_i in the natural basis
		tmpMatrix = findSxMatrices(0, natBasis);
		typename OperatorType::Su2RelatedType su2related3;
		OperatorType                          myOp3(tmpMatrix,
                                   ProgramGlobals::FermionOrBosonEnum::BOSON,
                                   PairType(2, 1),
                                   1.0 / sqrt(2.0),
                                   su2related3);
		this->createOpsLabel("sx").push(myOp3);
		this->makeTrackable("sx");
	}

	void fillModelLinks() override
	{
		if (BasisType::useSu2Symmetry())
			err("SU(2): no longer supported\n");

		// auto lambda = []
		ModelTermType& plaquetteX = ModelBaseType::createTerm("PlaquetteX");

		OpForLinkType sx("sx");

		plaquetteX.push4(sx, 'N', sx, 'N', sx, 'N', sx, 'N');
	}

	void fillNewNonLocals(std::vector<OperatorType>& newNonLocals,
	                      const LeftRightSuperType&  lrs,
	                      RealType) const override // last argument is time
	{
		bool isSys = (lrs.super().block()[0] == 0);
		if (isSys) {
			// (0, 1) (0, 1, 2) (2, 3), (2,3,4)...
			SizeType nsuper = lrs.super().block().size();
			assert(nsuper > 0);
			SizeType last = lrs.super().block()[nsuper - 1];
			assert(last > 0);
			assert(lrs.right().block().size() == 1);
			if (last & 1) {
				// 0 x 1 (or 2 x 3 etc.)
				const OperatorType& op1 = lrs.left().localOperator(last - 1);
				const OperatorType& op2 = lrs.right().localOperator(0);
				newNonLocals.push_back(op1 * op2);
			} else {
				// (0,1) x (2) (or (2,3)x(4) etc.)
				SizeType            lastNonLocal = last / 2 - 1;
				const OperatorType& op1 = lrs.left().getSuperByIndex(lastNonLocal);
				const OperatorType& op2 = lrs.right().localOperator(0);
				newNonLocals.push_back(op1 * op2);
			}
		} else {
			throw PsimagLite::RuntimeError("fillNewNonLocals unimplemented for env\n");
		}
	}

	SuperOpHelperBaseType* setSuperOpHelper() override
	{
		if (!superOpHelperPlaquette_)
			superOpHelperPlaquette_
			    = new SuperOpHelperPlaquetteType(ModelBaseType::superGeometry());
		return superOpHelperPlaquette_;
	}

private:

	//! Find S^+_site in the natural basis natBasis
	SparseMatrixType findSplusMatrices(SizeType site, const HilbertBasisType& natBasis) const
	{
		SizeType   total = natBasis.size();
		MatrixType cm(total, total);
		RealType   j              = 0.5 * TWICE_THE_SPIN;
		SizeType   bitsForOneSite = utils::bitSizeOfInteger(TWICE_THE_SPIN);
		SizeType   bits           = 1 + ProgramGlobals::logBase2(TWICE_THE_SPIN);
		SizeType   mask           = 1;
		mask <<= bits; // mask = 2^bits
		assert(mask > 0);
		mask--;
		mask <<= (site * bitsForOneSite);

		for (SizeType ii = 0; ii < total; ii++) {
			SizeType ket = natBasis[ii];

			SizeType ketsite = ket & mask;
			ketsite >>= (site * bitsForOneSite);

			assert(ketsite == ket);
			SizeType brasite = ketsite + 1;
			if (brasite >= TWICE_THE_SPIN + 1)
				continue;

			SizeType bra = ket & (~mask);
			assert(bra == 0);
			brasite <<= (site * bitsForOneSite);
			bra |= brasite;
			assert(bra == brasite);

			RealType m = ketsite - j;
			RealType x = j * (j + 1) - m * (m + 1);
			assert(x >= 0);

			cm(ket, bra) = sqrt(x);
		}

		SparseMatrixType operatorMatrix(cm);
		return operatorMatrix;
	}

	//! Find S^z_i in the natural basis natBasis
	SparseMatrixType findSzMatrices(SizeType site, const HilbertBasisType& natBasis) const
	{
		SizeType   total = natBasis.size();
		MatrixType cm(total, total);
		RealType   j              = 0.5 * TWICE_THE_SPIN;
		SizeType   bitsForOneSite = utils::bitSizeOfInteger(TWICE_THE_SPIN);
		SizeType   bits           = ProgramGlobals::logBase2(TWICE_THE_SPIN) + 1;
		SizeType   mask           = 1;
		mask <<= bits; // mask = 2^bits
		assert(mask > 0);
		mask--;
		mask <<= (site * bitsForOneSite);

		for (SizeType ii = 0; ii < total; ii++) {
			SizeType ket = natBasis[ii];

			SizeType ketsite = ket & mask;
			ketsite >>= (site * bitsForOneSite);
			assert(ketsite == ket);
			RealType m   = ketsite - j;
			cm(ket, ket) = m;
		}

		SparseMatrixType operatorMatrix(cm);
		return operatorMatrix;
	}

	SparseMatrixType findSxMatrices(SizeType site, const HilbertBasisType& natBasis) const
	{
		SparseMatrixType Splus_temp = findSplusMatrices(site, natBasis);
		SparseMatrixType Sminus_temp, Sx;
		transposeConjugate(Sminus_temp, Splus_temp);
		RealType tmp = 0.5;

		Sx = tmp * Splus_temp;
		Sx += tmp * Sminus_temp;

		return Sx;
	}

	void setSymmetryRelated(VectorQnType& qns, const HilbertBasisType& basis) const
	{
		// find j,m and flavors (do it by hand since we assume n==1)
		// note: we use 2j instead of j
		// note: we use m+j instead of m
		// This assures us that both j and m are SizeType
		using PairType = std::pair<SizeType, SizeType>;

		VectorSizeType other;
		QnType::ifPresentOther0IsElectrons = false;
		qns.resize(basis.size(), QnType::zero());
		for (SizeType i = 0; i < basis.size(); ++i) {
			PairType jmpair(0, 0);
			SizeType flavor = 1;
			qns[i]          = QnType(false, other, jmpair, flavor);
		}
	}

	GaugeSpin(const GaugeSpin&) = delete;

	GaugeSpin& operator=(const GaugeSpin&) = delete;

	ParametersGaugeSpin<RealType, QnType> modelParameters_;
	SuperOpHelperPlaquetteType*           superOpHelperPlaquette_;
}; // class GaugeSpin

} // namespace Dmrg
/*@}*/
#endif // DMRG_MODEL_GAUGESPIN_H
