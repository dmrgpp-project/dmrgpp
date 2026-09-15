// SPDX-FileCopyrightText: Copyright (c) 2009, 2017-2026, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 6+]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file LiouvillianHeisenberg.hh
 *
 *  An implementation of the LiouvilleHeisenberg Model to use with  DmrgSolver
 *
 */

#ifndef DMRG_LIOUVILLIANHEISENBERG_H
#define DMRG_LIOUVILLIANHEISENBERG_H

#include "ImaginaryUnitOrFail.hh"
#include "ParamsLiouvillianHeisenberg.hh"
#include "ProgramGlobals.h"
#include "Utils.h"
#include "VerySparseMatrix.h"
#include <PsimagLite/CrsMatrix.h>
#include <algorithm>
#include <utility>

namespace Dmrg {

template <typename ModelBaseType> class LiouvillianHeisenberg : public ModelBaseType {
public:

	using ModelHelperType      = typename ModelBaseType::ModelHelperType;
	using BasisType            = typename ModelHelperType::BasisType;
	using SuperGeometryType    = typename ModelBaseType::SuperGeometryType;
	using LeftRightSuperType   = typename ModelBaseType::LeftRightSuperType;
	using LinkType             = typename ModelBaseType::LinkType;
	using OperatorsType        = typename ModelHelperType::OperatorsType;
	using RealType             = typename ModelHelperType::RealType;
	using VectorType           = typename ModelBaseType::VectorType;
	using QnType               = typename ModelBaseType::QnType;
	using VectorQnType         = typename ModelBaseType::VectorQnType;
	using BlockType            = typename ModelBaseType::BlockType;
	using SolverParamsType     = typename ModelBaseType::SolverParamsType;
	using SparseMatrixType     = typename ModelHelperType::SparseMatrixType;
	using ComplexOrRealType    = typename SparseMatrixType::value_type;
	using WordType             = unsigned int long;
	using InputValidatorType   = typename ModelBaseType::InputValidatorType;
	using MatrixType           = PsimagLite::Matrix<ComplexOrRealType>;
	using VectorSizeType       = typename PsimagLite::Vector<SizeType>::Type;
	using VectorRealType       = typename ModelBaseType::VectorRealType;
	using ModelTermType        = typename ModelBaseType::ModelTermType;
	using HilbertBasisType     = typename PsimagLite::Vector<SizeType>::Type;
	using OperatorType         = typename OperatorsType::OperatorType;
	using PairType             = typename OperatorType::PairType;
	using VectorOperatorType   = typename PsimagLite::Vector<OperatorType>::Type;
	using MyBasis              = typename ModelBaseType::MyBasis;
	using MyBasisWithOperators = typename ModelBaseType::BasisWithOperatorsType;
	using OpsLabelType         = typename ModelBaseType::OpsLabelType;
	using OpForLinkType        = typename ModelBaseType::OpForLinkType;
	using ModelParametersType  = ParamsLiouvillianHeisenberg<RealType, QnType>;

	LiouvillianHeisenberg(const SolverParamsType&  solverParams,
	                      InputValidatorType&      io,
	                      const SuperGeometryType& geometry)
	    : ModelBaseType(solverParams, geometry, io)
	    , modelParameters_(io)
	    , superGeometry_(geometry)
	    , is_separate_(modelParameters_.ancillas == ModelParametersType::Ancillas::SEPARATE)

	{
		SizeType n = superGeometry_.numberOfSites();

		ModelParametersType::checkMagneticField(
		    modelParameters_.magneticFieldZ.size(), 'Z', n);
	}

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

		assert(block.size() == 1);
		SizeType site = block[0];
		addMagneticFieldZ(hmatrix, modelParameters_.magneticFieldZ, site);

		// physical ancilla jump connection
		const auto iter = jump_mapping_.find(site);
		if (iter == jump_mapping_.end()) {
			// ATTENTION: Early exit here
			return;
		}

		SizeType index = iter->second;
		assert(index < jump_operators_.size());

		hmatrix += jump_operators_[index];
	}

	bool isHermitian() const final { return false; }

protected:

	void fillLabeledOperators(VectorQnType& qns) override
	{
		SizeType hilbert_small = modelParameters_.twiceTheSpin + 1; // = 2
		SizeType hilbert_big
		    = hilbert_small * hilbert_small; // physical and ancilla together
		SizeType total = (is_separate_) ? hilbert_small : hilbert_big;

		std::vector<SizeType> perm(total);
		{
			// Do not use natBasis outside this block
			HilbertBasisType natBasis(total);
			for (SizeType i = 0; i < total; ++i) {
				natBasis[i] = i;
				perm[i]     = i;
			}

			setSymmetryRelated(qns, natBasis);
		}

		std::vector<RealType> signs(hilbert_small, 1);

		// Small size operators
		SparseMatrixType splus_small = findSplusMatrices(); // 2-state space only
		SparseMatrixType sz_small    = findSzMatrices(); // 2-state space only

		// Set the operators S^+ in the natural basis
		SparseMatrixType tmpMatrix;
		std::string      label;

		if (is_separate_) {
			label     = "splus";
			tmpMatrix = splus_small;
		} else {
			// Physical S+
			label = "splus_p";
			externalProduct(tmpMatrix, splus_small, hilbert_small, signs, false, perm);
		}

		tmpMatrix.checkValidity();
		OperatorType myOp1(tmpMatrix, ProgramGlobals::FermionOrBosonEnum::BOSON);
		this->createOpsLabel(label).push(myOp1);
		this->makeTrackable(label);

		// Physical S-
		myOp1.dagger();
		label = (is_separate_) ? "sminus" : "sminus_p";
		this->createOpsLabel(label).push(myOp1);

		if (!is_separate_) {
			// Ancilla S+ and S-
			externalProduct(tmpMatrix, splus_small, hilbert_small, signs, true, perm);
			tmpMatrix.checkValidity();
			OperatorType myOp2(tmpMatrix, ProgramGlobals::FermionOrBosonEnum::BOSON);
			this->createOpsLabel("splus_a").push(myOp2);
			this->makeTrackable("splus_a");

			myOp2.dagger();
			this->createOpsLabel("sminus_a").push(myOp2);
		}

		// Set the operators S^z in the natural basis
		if (is_separate_) {
			label     = "sz";
			tmpMatrix = sz_small;
		} else {
			// Physical Sz
			label = "sz_p";
			externalProduct(tmpMatrix, sz_small, hilbert_small, signs, false, perm);
		}

		tmpMatrix.checkValidity();
		OperatorType myOp3(tmpMatrix, ProgramGlobals::FermionOrBosonEnum::BOSON);
		this->createOpsLabel(label).push(myOp3);
		this->makeTrackable(label);

		if (!is_separate_) {
			// Ancilla Sz
			externalProduct(tmpMatrix, sz_small, hilbert_small, signs, true, perm);
			tmpMatrix.checkValidity();
			OperatorType myOp4(tmpMatrix, ProgramGlobals::FermionOrBosonEnum::BOSON);
			this->createOpsLabel("sz_a").push(myOp4);
			this->makeTrackable("sz_a");
		}

		cacheJumpOperators();
	}

	void fillModelLinks() override
	{
		// physical H
		// Multiply by -i
		ComplexOrRealType sqrt_minus_one = ImginaryUnitOrFail<ComplexOrRealType>::value();
		std::string       addition       = (is_separate_) ? "" : "_p";
		connectionsSpSm("spsm" + addition, "splus" + addition, -sqrt_minus_one);
		connectionSzSz("szsz" + addition, "sz" + addition, -sqrt_minus_one);

		if (!is_separate_) {
			// ancilla H^T
			// Multiply by +i
			connectionsSpSm("spsm_a", "splus_a", sqrt_minus_one);
			connectionSzSz("szsz_a", "sz_a", sqrt_minus_one);
		} else {
			// jump operators become connections when phys and anc are separate

			ModelTermType& spsp = ModelBaseType::createTerm("spsp");
			OpForLinkType  splus("splus");
			spsp.push(splus, 'N', splus, 'N');

			ModelTermType& smsm = ModelBaseType::createTerm("smsm");
			smsm.push(splus, 'C', splus, 'C');
		}
	}

private:

	void addMagneticFieldZ(SparseMatrixType&            hmatrix,
	                       const std::vector<RealType>& v,
	                       SizeType                     site) const
	{
		const SizeType linSize = superGeometry_.numberOfSites();
		if (v.size() != linSize) {
			// no magnetic field supplied
			if (v.size() == 0)
				return;

			// magnetic field supplied with wrong size
			// should have been caught earlier
			assert(false);
		}

		assert(site < v.size());
		RealType            tmp   = v[site];
		std::string         label = (is_separate_) ? "sz" : "sz_p";
		const OperatorType& sz_p  = ModelBaseType::naturalOperator(label, site, 0);
		// MagneticFieldZ is the coefficient of Sz in H, not in L.
		// For separate sites, supply signed entries: +b on the physical
		// leg and -b on its ancilla, as for the exchange connectors.
		const ComplexOrRealType imaginaryUnit
		    = ImginaryUnitOrFail<ComplexOrRealType>::value();
		hmatrix += (-imaginaryUnit * tmp) * sz_p.getCRS();

		if (!is_separate_) {
			const OperatorType& sz_a = ModelBaseType::naturalOperator("sz_a", site, 0);
			hmatrix += (imaginaryUnit * tmp) * sz_a.getCRS();
		}
	}

	void connectionsSpSm(const std::string&       connection_name,
	                     const std::string&       op_name,
	                     const ComplexOrRealType& factor)
	{
		// The commutator prefactor (+/- i) is the same for both exchange
		// directions. Automatic Hermitian completion would conjugate it.
		ModelTermType& spsm = ModelBaseType::createTerm(connection_name, false);
		OpForLinkType  splus(op_name);

		auto valueModiferTerm0
		    = [factor](ComplexOrRealType& value) { value *= (0.5 * factor); };
		spsm.push(splus, 'N', splus, 'C', valueModiferTerm0);
		spsm.push(splus, 'C', splus, 'N', valueModiferTerm0);
	}

	void connectionSzSz(const std::string&       connection_name,
	                    const std::string&       op_name,
	                    const ComplexOrRealType& factor)
	{
		ModelTermType& szsz = ModelBaseType::createTerm(connection_name);

		OpForLinkType sz(op_name);
		auto valueModiferTerm0 = [factor](ComplexOrRealType& value) { value *= factor; };
		szsz.push(sz, 'N', sz, 'N', valueModiferTerm0);
	}

	//! Find S^+_site in the physical "small" basis only
	SparseMatrixType findSplusMatrices() const
	{
		// We assume there's only one kind of sites
		// FIXME TODO for SDHS
		// SDHS = site-dependent Hilbert spaces
		SizeType site = 0;

		SizeType total = modelParameters_.twiceTheSpin + 1;
		assert(total == 2);
		MatrixType cm(total, total);
		RealType   j              = 0.5 * modelParameters_.twiceTheSpin;
		SizeType   bitsForOneSite = utils::bitSizeOfInteger(modelParameters_.twiceTheSpin);
		SizeType   bits = 1 + ProgramGlobals::logBase2(modelParameters_.twiceTheSpin);
		SizeType   mask = 1;
		mask <<= bits; // mask = 2^bits
		assert(mask > 0);
		mask--;
		mask <<= (site * bitsForOneSite);

		for (SizeType ii = 0; ii < total; ii++) {
			SizeType ket = ii;

			SizeType ketsite = ket & mask;
			ketsite >>= (site * bitsForOneSite);

			assert(ketsite == ket);
			SizeType brasite = ketsite + 1;
			if (brasite >= modelParameters_.twiceTheSpin + 1)
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

	//! Find S^z_i in the physical small basis only
	SparseMatrixType findSzMatrices() const
	{
		// We assume there's only one kind of sites
		// FIXME TODO for SDHS
		// SDHS = site-dependent Hilbert spaces
		SizeType site = 0;

		SizeType total = modelParameters_.twiceTheSpin + 1;
		assert(total == 2);
		MatrixType cm(total, total);
		RealType   j              = 0.5 * modelParameters_.twiceTheSpin;
		SizeType   bitsForOneSite = utils::bitSizeOfInteger(modelParameters_.twiceTheSpin);
		SizeType   bits = ProgramGlobals::logBase2(modelParameters_.twiceTheSpin) + 1;
		SizeType   mask = 1;
		mask <<= bits; // mask = 2^bits
		assert(mask > 0);
		mask--;
		mask <<= (site * bitsForOneSite);

		for (SizeType ii = 0; ii < total; ii++) {
			SizeType ket = ii;

			SizeType ketsite = ket & mask;
			ketsite >>= (site * bitsForOneSite);
			assert(ketsite == ket);
			RealType m   = ketsite - j;
			cm(ket, ket) = m;
		}

		SparseMatrixType operatorMatrix(cm);
		return operatorMatrix;
	}

	// There are no symmetries due to the s^- x s^- term and similar for s^+ s^+
	void setSymmetryRelated(VectorQnType& qns, const HilbertBasisType& basis) const
	{
		using PairType = std::pair<SizeType, SizeType>;

		VectorSizeType other(1);
		QnType::ifPresentOther0IsElectrons = false;
		qns.resize(basis.size(), QnType::zero());
		for (SizeType i = 0; i < basis.size(); ++i) {
			SizeType flavor = 1;
			qns[i]          = QnType(false, other, PairType(0, 0), flavor);
		}
	}

	void cacheJumpOperators()
	{
		SizeType n = superGeometry_.numberOfSites();
		if (modelParameters_.bath_gamma.size() != n) {
			err("bath_gamma.size need to be equal to number of sites\n");
		}

		if (modelParameters_.bath_f.size() != n) {
			err("bath_f.size need to be equal to number of sites\n");
		}

		for (SizeType site = 0; site < n; ++site) {
			if (modelParameters_.bath_gamma[site] == 0)
				continue;

			ComplexOrRealType factor_a
			    = modelParameters_.bath_gamma[site] * modelParameters_.bath_f[site];
			ComplexOrRealType factor_b = modelParameters_.bath_gamma[site] - factor_a;

			MatrixType dense = (is_separate_)
			    ? jumpOperatorWhenSeparate(site, factor_a, factor_b)
			    : jumpOperatorWhenTogether(site, factor_a, factor_b);

			SparseMatrixType tmpMatrix(dense);
			jump_mapping_[site] = jump_operators_.size();
			jump_operators_.push_back(tmpMatrix);
		}
	}

	// When separate, in the input file, if a physical site has f_i its ancilla ought to
	// have 1 - f_i
	MatrixType jumpOperatorWhenSeparate(SizeType                 site,
	                                    const ComplexOrRealType& factor_a,
	                                    const ComplexOrRealType& factor_b) const
	{
		SparseMatrixType sm = ModelBaseType::naturalOperator("sminus", site, 0).getCRS();
		SparseMatrixType sp = ModelBaseType::naturalOperator("splus", site, 0).getCRS();

		SparseMatrixType p0 = sm * sp;
		SparseMatrixType p1 = sp * sm;

		MatrixType p0_dense = p0.toDense();
		MatrixType p1_dense = p1.toDense();

		SizeType rank = p0_dense.rows();
		assert(rank == p0_dense.cols());
		assert(rank == 2);
		MatrixType tmp_dense(rank, rank);
		for (SizeType i = 0; i < rank; ++i) {
			for (SizeType j = 0; j < rank; ++j) {
				tmp_dense(i, j) = -0.5 * factor_a * p1_dense(i, j)
				    - 0.5 * factor_b * p0_dense(i, j);
			}
		}

		return tmp_dense;
	}

	MatrixType jumpOperatorWhenTogether(SizeType                 site,
	                                    const ComplexOrRealType& factor_a,
	                                    const ComplexOrRealType& factor_b) const
	{
		SparseMatrixType sm_p
		    = ModelBaseType::naturalOperator("sminus_p", site, 0).getCRS();
		SparseMatrixType sm_a
		    = ModelBaseType::naturalOperator("sminus_a", site, 0).getCRS();
		SparseMatrixType sp_p = ModelBaseType::naturalOperator("splus_p", site, 0).getCRS();
		SparseMatrixType sp_a = ModelBaseType::naturalOperator("splus_a", site, 0).getCRS();

		SparseMatrixType p0_p = sm_p * sp_p;
		SparseMatrixType p0_a = sm_a * sp_a;
		SparseMatrixType p1_p = sp_p * sm_p;
		SparseMatrixType p1_a = sp_a * sm_a;

		SizeType hilbert = sm_p.rows();
		assert(hilbert == 4); // because of the ancillas
		SizeType physical_hilbert = sqrt(hilbert);
		assert(physical_hilbert == 2); // spin 1/2
		std::vector<RealType> signs(physical_hilbert, 1);
		std::vector<SizeType> perm(hilbert);
		for (SizeType i = 0; i < hilbert; ++i) {
			perm[i] = i;
		}

		SparseMatrixType sp = findSplusMatrices();
		SparseMatrixType sm;
		transposeConjugate(sm, sp);

		SparseMatrixType sm_cross_sm;
		SparseMatrixType sp_cross_sp;
		externalProduct(sm_cross_sm, sm, sm, signs, true, perm);
		externalProduct(sp_cross_sp, sp, sp, signs, true, perm);

		MatrixType       sm_cross_sm_dense = sm_cross_sm.toDense();
		MatrixType       sp_cross_sp_dense = sp_cross_sp.toDense();
		SparseMatrixType p0                = p0_p;
		p0 += p0_a;
		SparseMatrixType p1 = p1_p;
		p1 += p1_a;
		MatrixType p0_dense = p0.toDense();
		MatrixType p1_dense = p1.toDense();

		SizeType rank = p0_dense.rows();
		assert(rank == p0_dense.cols());
		assert(rank == 4);
		MatrixType tmp_dense(rank, rank);
		for (SizeType i = 0; i < rank; ++i) {
			for (SizeType j = 0; j < rank; ++j) {
				tmp_dense(i, j)
				    = factor_a * (sm_cross_sm_dense(i, j) - 0.5 * p1_dense(i, j))
				    + factor_b * (sp_cross_sp_dense(i, j) - 0.5 * p0_dense(i, j));
			}
		}

		return tmp_dense;
	}

	ModelParametersType           modelParameters_;
	const SuperGeometryType&      superGeometry_;
	bool                          is_separate_;
	std::vector<SparseMatrixType> jump_operators_;
	std::map<SizeType, SizeType>  jump_mapping_;
}; // class LiouvillianHeisenberg

} // namespace Dmrg
/*@}*/
#endif // DMRG_LIOUVILLIANHEISENBERG_H
