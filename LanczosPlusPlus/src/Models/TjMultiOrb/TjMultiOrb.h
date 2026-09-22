/*
 */

#ifndef LANCZOS_TJ_MULTIORB_H
#define LANCZOS_TJ_MULTIORB_H

#include "BasisTjMultiOrbLanczos.h"
#include "LanczosModelBase.hpp"
#include "ParametersTjMultiOrb.h"
#include <PsimagLite/BitManip.h>
#include <PsimagLite/CrsMatrix.h>
#include <PsimagLite/SparseRow.h>
#include <PsimagLite/TypeToString.h>

namespace LanczosPlusPlus {

template <typename ComplexOrRealType, typename GeometryType, typename InputType>
class TjMultiOrb : public LanczosModelBase<ComplexOrRealType, GeometryType, InputType> {

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type           RealType;
	typedef PsimagLite::Matrix<ComplexOrRealType>                        MatrixType;
	typedef std::pair<SizeType, SizeType>                                PairType;
	typedef LanczosModelBase<ComplexOrRealType, GeometryType, InputType> BaseType;

	enum
	{
		SPIN_UP   = LanczosGlobals::SPIN_UP,
		SPIN_DOWN = LanczosGlobals::SPIN_DOWN
	};

	static const SizeType REINTERPRET_6 = 6;
	static const SizeType REINTERPRET_9 = 9;
	static const SizeType STATE_EMPTY   = 0;
	static const SizeType STATE_UP_A    = 1;
	static const SizeType STATE_DOWN_A  = 4;

public:

	typedef ParametersTjMultiOrb<RealType, InputType>       ParametersModelType;
	typedef BasisTjMultiOrbLanczos<GeometryType>            BasisType;
	typedef typename BasisType::BaseType                    BasisBaseType;
	typedef typename BasisType::WordType                    WordType;
	typedef typename BaseType::SparseMatrixType             SparseMatrixType;
	typedef typename BaseType::VectorType                   VectorType;
	typedef typename BaseType::VectorSizeType               VectorSizeType;
	typedef std::pair<WordType, WordType>                   PairWordType;
	typedef typename PsimagLite::Vector<PairWordType>::Type VectorPairWordType;
	typedef PsimagLite::Matrix<SizeType>                    MatrixSizeType;
	typedef PsimagLite::SparseRow<SparseMatrixType>         SparseRowType;

	TjMultiOrb(SizeType nup, SizeType ndown, InputType& io, const GeometryType& geometry)
	    : mp_(io)
	    , geometry_(geometry)
	    , basis_(geometry, nup, ndown, mp_.orbitals)
	    , hoppings_(geometry_.numberOfSites() * geometry_.numberOfSites(),
	                mp_.orbitals * mp_.orbitals)
	    , jpm_(geometry_.numberOfSites() * geometry_.numberOfSites(),
	           mp_.orbitals * mp_.orbitals)
	    , jzz_(geometry_.numberOfSites() * geometry_.numberOfSites(),
	           mp_.orbitals * mp_.orbitals)
	    , w_(geometry_.numberOfSites() * geometry_.numberOfSites(), mp_.orbitals * mp_.orbitals)
	{
		SizeType n = geometry_.numberOfSites();

		if (geometry_.terms() != 4) {
			PsimagLite::String msg("TjMultiOrb: must have 4 terms\n");
			throw PsimagLite::RuntimeError(msg);
		}

		for (SizeType i = 0; i < n; i++) {
			for (SizeType j = 0; j < n; j++) {
				for (SizeType orb1 = 0; orb1 < mp_.orbitals; ++orb1) {
					for (SizeType orb2 = 0; orb2 < mp_.orbitals; ++orb2) {
						hoppings_(i + j * n, orb1 + orb2 * mp_.orbitals)
						    = geometry_(i, orb1, j, orb2, 0);
						jpm_(i + j * n, orb1 + orb2 * mp_.orbitals)
						    = geometry_(i, orb1, j, orb2, 1);
						jzz_(i + j * n, orb1 + orb2 * mp_.orbitals)
						    = geometry_(i, orb1, j, orb2, 2);
						w_(i + j * n, orb1 + orb2 * mp_.orbitals)
						    = geometry_(i, orb1, j, orb2, 3);
					}
				}
			}
		}
	}

	~TjMultiOrb() override { BaseType::deleteGarbage(garbage_); }

	SizeType size() const override { return basis_.size(); }

	SizeType orbitals(SizeType) const override { return mp_.orbitals; }

	void setupHamiltonian(SparseMatrixType& matrix) const override
	{
		setupHamiltonian(matrix, basis_);
	}

	//! Gf. related functions below:
	void setupHamiltonian(SparseMatrixType& matrix, const BasisBaseType& basis) const override
	{
		SizeType                                    hilbert = basis.size();
		typename PsimagLite::Vector<RealType>::Type diag(hilbert, 0.0);
		calcDiagonalElements(diag, basis);

		SizeType nsite = geometry_.numberOfSites();

		matrix.resize(hilbert, hilbert);
		// Calculate off-diagonal elements AND store matrix
		SizeType n_counter = 0;
		for (SizeType ispace = 0; ispace < hilbert; ispace++) {
			SparseRowType sparse_row;
			matrix.setRow(ispace, n_counter);
			WordType ket1 = basis(ispace, SPIN_UP);
			WordType ket2 = basis(ispace, SPIN_DOWN);
			// Save diagonal
			sparse_row.add(ispace, diag[ispace]);
			for (SizeType i = 0; i < nsite; i++) {
				for (SizeType orb = 0; orb < mp_.orbitals; ++orb) {
					setHoppingTerm(sparse_row, ket1, ket2, i, orb, basis);
				}
				setSplusSminus(sparse_row, ket1, ket2, i, basis);
			}

			n_counter += sparse_row.finalize(matrix);
		}

		matrix.setRow(hilbert, n_counter);
		matrix.checkValidity();
		assert(isHermitian(matrix));
		//		std::cout<<"MATRIX before rotation\n";
		//		std::cout<<matrix.toDense();
		if (mp_.reinterpretAndTruncate)
			reinterpretAndTruncate(matrix, basis);
		//		std::cout<<"MATRIX after rotation and truncation\n";
		//		std::cout<<matrix.toDense();
	}

	bool hasNewParts(std::pair<SizeType, SizeType>&       newParts,
	                 const std::pair<SizeType, SizeType>& oldParts,
	                 const LabeledOperator&               lOperator,
	                 SizeType                             spin,
	                 SizeType) const override
	{
		if (lOperator.id() == LabeledOperator::Label::OPERATOR_C
		    || lOperator.id() == LabeledOperator::Label::OPERATOR_CDAGGER)
			return hasNewPartsCorCdagger(newParts, oldParts, lOperator, spin);

		if (lOperator.id() == LabeledOperator::Label::OPERATOR_SPLUS
		    || lOperator.id() == LabeledOperator::Label::OPERATOR_SMINUS)
			return hasNewPartsSplusOrMinus(newParts, oldParts, lOperator, spin);

		PsimagLite::String str(__FILE__);
		str += " " + ttos(__LINE__) + "\n";
		str += PsimagLite::String("hasNewParts: unsupported operator ");
		str += lOperator.toString() + "\n";
		throw std::runtime_error(str.c_str());
	}

	const GeometryType& geometry() const override { return geometry_; }

	const BasisType& basis() const override { return basis_; }

	void printBasis(std::ostream& os) const { os << basis_; }

	PsimagLite::String name() const override { return __FILE__; }

	BasisType* createBasis(SizeType nup, SizeType ndown) const override
	{
		BasisType* ptr = new BasisType(geometry_, nup, ndown, mp_.orbitals);
		garbage_.push_back(ptr);
		return ptr;
	}

	void
	matrixVectorProduct(VectorType&, const VectorType&, const BasisBaseType&) const override
	{
		throw PsimagLite::RuntimeError(
		    "ModelBase::matrixVectorProduct(3) not impl. for this model\n");
	}

	void print(std::ostream& os) const override { os << mp_; }

	void printOperators(std::ostream& os) const override
	{
		SizeType nup   = basis_.electrons(SPIN_UP);
		SizeType ndown = basis_.electrons(SPIN_DOWN);
		os << "#SectorSource 2 " << nup << " " << ndown << "\n";
		SizeType spin = SPIN_UP;
		for (SizeType site = 0; site < geometry_.numberOfSites(); ++site)
			printOperatorC(site, spin, os);
	}

private:

	void reinterpretAndTruncate(SparseMatrixType& matrix, const BasisBaseType& basis) const
	{
		SizeType         n = basis.size();
		SparseMatrixType rot(n, n);
		SparseMatrixType rot_t;
		VectorSizeType   targets;
		buildRotation(rot, rot_t, targets, basis);
		//		std::cout<<"Rotation\n";
		//		std::cout<<rot.toDense();
		SparseMatrixType tmp;
		multiply(tmp, matrix, rot_t);
		multiply(matrix, rot, tmp);
		assert(isHermitian(matrix));
		//		std::cout<<"Matrix after rotation\n";
		//		std::cout<<matrix.toDense();
		truncateMatrix(matrix, targets);
		assert(isHermitian(matrix));
	}

	void truncateMatrix(SparseMatrixType& matrix, const VectorSizeType& targets) const
	{
		assert(matrix.rows() > targets.size());
		SizeType         n_full  = matrix.rows();
		SizeType         n_trunc = matrix.rows() - targets.size();
		SparseMatrixType matrix2(n_trunc, n_trunc);

		VectorSizeType remap(n_full, 0);
		SizeType       ii = 0;
		for (SizeType i = 0; i < n_full; ++i) {
			if (std::find(targets.begin(), targets.end(), i) != targets.end())
				continue;
			remap[i] = ii;
			ii++;
		}

		assert(ii == nTrunc);
		SizeType counter = 0;
		for (SizeType i = 0; i < n_full; ++i) {
			if (std::find(targets.begin(), targets.end(), i) != targets.end())
				continue;
			matrix2.setRow(remap[i], counter);
			for (int k = matrix.getRowPtr(i); k < matrix.getRowPtr(i + 1); ++k) {
				SizeType col = matrix.getCol(k);
				if (std::find(targets.begin(), targets.end(), col) != targets.end())
					continue;
				matrix2.pushCol(remap[col]);
				matrix2.pushValue(matrix.getValue(k));
				counter++;
			}
		}

		matrix2.setRow(n_trunc, counter);
		matrix2.checkValidity();
		matrix = matrix2;
	}

	void buildRotation(SparseMatrixType&    rot,
	                   SparseMatrixType&    rotT,
	                   VectorSizeType&      targets,
	                   const BasisBaseType& basis) const
	{
		SizeType n_counter = 0;
		SizeType hilbert   = basis.size();
		for (SizeType ispace = 0; ispace < hilbert; ispace++) {
			SparseRowType sparse_row;
			rot.setRow(ispace, n_counter);
			WordType       ket1 = basis(ispace, SPIN_UP);
			WordType       ket2 = basis(ispace, SPIN_DOWN);
			VectorSizeType k;
			breakIntoSites(k, ket1, ket2);
			SizeType       branches = getBranchesFromVector(k);
			VectorType     bra_values(branches, 1.0);
			MatrixSizeType bra_matrix(branches, k.size());
			getBrasFromVector(bra_matrix, bra_values, k);
			for (SizeType b = 0; b < branches; ++b) {
				PairWordType bra  = fromMatrixToBra(bra_matrix, b);
				SizeType     temp = basis.perfectIndex(bra.first, bra.second);
				sparse_row.add(temp, bra_values[b]);
				if (isToBeRemoved(bra.first, bra.second))
					targets.push_back(temp);
			}

			n_counter += sparse_row.finalize(rot);
		}

		rot.setRow(hilbert, n_counter);
		rot.checkValidity();
		transposeConjugate(rotT, rot);

		PsimagLite::Sort<VectorSizeType> sort;
		VectorSizeType                   iperm(targets.size(), 0);
		sort.sort(targets, iperm);
		typename VectorSizeType::iterator it = std::unique(targets.begin(), targets.end());
		SizeType                          new_size = it - targets.begin();
		targets.resize(new_size);
	}

	bool isToBeRemoved(WordType bra1, WordType bra2) const
	{
		SizeType       n = geometry_.numberOfSites();
		VectorSizeType k(n, 0);
		breakIntoSitesEx(k, bra1, bra2);
		assert(k.size() == geometry_.numberOfSites());
		for (SizeType i = 0; i < k.size(); ++i)
			if (isToBeRemoved(k[i]))
				return true;

		return false;
	}

	bool isToBeRemoved(SizeType one) const
	{
		assert(mp_.reinterpretAndTruncate > 0);
		bool b = false;
		if (mp_.reinterpretAndTruncate > 0)
			b |= hasState(one, REINTERPRET_6);
		if (mp_.reinterpretAndTruncate > 1)
			b |= hasState(one, STATE_EMPTY);
		if (mp_.reinterpretAndTruncate > 2) {
			b |= hasState(one, STATE_UP_A);
			b |= hasState(one, STATE_DOWN_A);
		}

		return b;
	}

	bool hasState(SizeType one, SizeType state) const { return (one == state); }

	SizeType getBranchesFromVector(const VectorSizeType& k) const
	{
		SizeType res = 1;
		for (SizeType i = 0; i < k.size(); ++i)
			if (k[i] == REINTERPRET_6 || k[i] == REINTERPRET_9)
				res *= 2;

		return res;
	}

	void breakIntoSites(VectorSizeType& k, WordType ket1, WordType ket2) const
	{
		SizeType mask = 3;
		while (ket1 > 0 || ket2 > 0) {
			SizeType val1 = ket1 & mask;
			SizeType val2 = ket2 & mask;
			val2 <<= 2;
			val1 |= val2;
			k.push_back(val1);
			ket1 >>= 2;
			ket2 >>= 2;
		}
	}

	void breakIntoSitesEx(VectorSizeType& k, WordType ket1, WordType ket2) const
	{
		SizeType n = k.size();
		assert(n == geometry_.numberOfSites());
		for (SizeType i = 0; i < n; ++i) {
			WordType one1 = ket1 & 3;
			WordType one2 = ket2 & 3;
			one2 <<= 2;
			k[i] = (one1 | one2);
			ket1 >>= 2;
			ket2 >>= 2;
			if (ket1 == 0 && ket2 == 0)
				break;
		}
	}
	// Rotation(old basis state) = linear combination of old basis states
	// The linear combination of old basis states contains various terms
	// These terms are called branches: the index b below
	// For each term we have a sites index
	// braMatrix: rows are terms or branch index, and cols are sites
	void getBrasFromVector(MatrixSizeType&       braMatrix,
	                       VectorType&           braValues,
	                       const VectorSizeType& k) const
	{
		SizeType   n                = k.size();
		SizeType   b                = 0;
		SizeType   current_branches = 1;
		SizeType   j                = 0;
		RealType   one_over_sqrt2   = 1.0 / sqrt(2.0);
		MatrixType m_values(braMatrix.n_row(), braMatrix.n_col());
		m_values.setTo(1.0);

		for (SizeType i = 0; i < n; ++i) {
			switch (k[i]) {
			// we're expanding old 6  = (new 6 + new 9)/sqrt(2)
			case REINTERPRET_6:
				for (b = 0; b < current_branches; ++b) {
					braMatrix(b, i)                    = REINTERPRET_6;
					braMatrix(current_branches + b, i) = REINTERPRET_9;
					m_values(b, i) *= one_over_sqrt2;
					m_values(current_branches + b, i) *= one_over_sqrt2;

					for (j = 0; j < i; ++j) {
						braMatrix(current_branches + b, j)
						    = braMatrix(b, j);
						m_values(current_branches + b, j) = m_values(b, j);
					}
				}

				current_branches *= 2;
				break;
				// we're expanding old 9  = (new 6 - new 9)/sqrt(2)
			case REINTERPRET_9:
				for (b = 0; b < current_branches; ++b) {
					braMatrix(b, i)                    = REINTERPRET_6;
					braMatrix(current_branches + b, i) = REINTERPRET_9;
					m_values(b, i) *= one_over_sqrt2;
					m_values(current_branches + b, i) *= (-one_over_sqrt2);

					for (j = 0; j < i; ++j) {
						braMatrix(current_branches + b, j)
						    = braMatrix(b, j);
						m_values(current_branches + b, j) = m_values(b, j);
					}
				}

				current_branches *= 2;
				break;
			default:
				for (b = 0; b < current_branches; ++b) {
					braMatrix(b, i) = k[i];
					m_values(b, i)  = 1.0;
				}

				break;
			}
		}

		assert(currentBranches == braMatrix.n_row());

		for (SizeType b = 0; b < m_values.n_row(); ++b) {
			ComplexOrRealType prod = 1;
			for (SizeType i = 0; i < m_values.n_col(); ++i)
				prod *= m_values(b, i);
			braValues[b] = prod;
		}
	}

	PairWordType fromMatrixToBra(const MatrixSizeType& braMatrix, SizeType branch) const
	{
		SizeType bra1  = 0;
		SizeType bra2  = 0;
		SizeType mask1 = 3;
		SizeType mask2 = 12;
		for (SizeType i = 0; i < braMatrix.n_col(); ++i) {
			SizeType tmp  = braMatrix(branch, i);
			SizeType tmp1 = tmp & mask1;
			tmp1 <<= (i * 2);
			bra1 |= tmp1;

			SizeType tmp2 = tmp & mask2;
			tmp2 >>= 2;
			tmp2 <<= (i * 2);
			bra2 |= tmp2;
		}

		return PairWordType(bra1, bra2);
	}

	void printOperatorC(SizeType site, SizeType spin, std::ostream& os) const
	{
		SizeType nup   = basis_.electrons(SPIN_UP);
		SizeType ndown = basis_.electrons(SPIN_DOWN);
		if (nup == 0) {
			os << "#Operator_c_" << spin << "_" << site << "\n";
			os << "#SectorDest 0\n"; // bogus
			os << "#Matrix\n";
			os << "0 0\n";
			return;
		}

		BasisType*     basis = createBasis(nup - 1, ndown);
		VectorSizeType opt(2, 0);
		opt[0] = site;
		opt[1] = spin;
		MatrixType matrix;
		setupOperator(matrix, *basis, "c", opt);
		os << "#Operator_c_" << spin << "_" << site << "\n";
		os << "#SectorDest 2 " << (nup - 1) << " " << ndown << "\n";
		os << "#Matrix\n";
		os << matrix;
	}

	void setupOperator(MatrixType&           matrix,
	                   const BasisBaseType&  basis,
	                   PsimagLite::String    operatorName,
	                   const VectorSizeType& operatorOptions) const
	{
		if (mp_.orbitals != 1)
			return;
		SizeType        hilbert_dest = basis.size();
		SizeType        hilbert_src  = basis_.size();
		SizeType        nsite        = geometry_.numberOfSites();
		LabeledOperator l_operator(operatorName);
		if (operatorName != "c") {
			PsimagLite::String str(__FILE__);
			str += " " + ttos(__LINE__) + "\n";
			str += "operator " + operatorName + " is unimplemented for this model\n";
			throw PsimagLite::RuntimeError(str);
		}

		if (operatorOptions.size() < 2) {
			PsimagLite::String str(__FILE__);
			str += " " + ttos(__LINE__) + "\n";
			str += "operator " + operatorName + " needs at least two options\n";
			throw PsimagLite::RuntimeError(str);
		}

		SizeType site = operatorOptions[0];
		if (site >= nsite) {
			PsimagLite::String str(__FILE__);
			str += " " + ttos(__LINE__) + "\n";
			str += "site requested " + ttos(site);
			str += " but number of sites= " + ttos(nsite) + "\n";
			throw PsimagLite::RuntimeError(str);
		}

		SizeType spin = operatorOptions[1];
		matrix.resize(hilbert_src, hilbert_dest);
		matrix.setTo(0.0);
		SizeType orb = 0;

		for (SizeType ispace = 0; ispace < hilbert_src; ispace++) {
			WordType ket1 = basis_(ispace, SPIN_UP);
			WordType ket2 = basis_(ispace, SPIN_DOWN);
			WordType bra  = ket1;
			// assumes OPERATOR_C
			bool b = basis.getBra(bra, ket1, ket2, l_operator, site, spin);
			if (!b)
				continue;
			SizeType index = basis.perfectIndex(bra, ket2);

			matrix(ispace, index) = basis.doSignGf(bra, ket2, site, spin, orb);
		}
	}

	bool hasNewPartsCorCdagger(std::pair<SizeType, SizeType>&       newParts,
	                           const std::pair<SizeType, SizeType>& oldParts,
	                           const LabeledOperator&               lOperator,
	                           SizeType                             spin) const
	{
		int new_part1 = oldParts.first;
		int new_part2 = oldParts.second;
		int c = (lOperator.id() == LabeledOperator::Label::OPERATOR_CDAGGER) ? 1 : -1;
		if (spin == SPIN_UP)
			new_part1 += c;
		else
			new_part2 += c;

		if (new_part1 < 0 || new_part2 < 0)
			return false;
		SizeType nsite = geometry_.numberOfSites();
		if (SizeType(new_part1) > nsite || SizeType(new_part2) > nsite)
			return false;
		if (new_part1 == 0 && new_part2 == 0)
			return false;
		if (SizeType(new_part1 + new_part2) > nsite)
			return false; // no double occupancy
		newParts.first  = SizeType(new_part1);
		newParts.second = SizeType(new_part2);
		return true;
	}

	bool hasNewPartsSplusOrMinus(std::pair<SizeType, SizeType>&       newParts,
	                             const std::pair<SizeType, SizeType>& oldParts,
	                             const LabeledOperator&               lOperator,
	                             SizeType                             spin) const
	{
		int new_part1 = oldParts.first;
		int new_part2 = oldParts.second;
		int c         = (lOperator.id() == LabeledOperator::Label::OPERATOR_SPLUS) ? 1 : -1;
		if (spin == SPIN_UP) {
			new_part1 += c;
			new_part2 -= c;
		} else {
			new_part2 += c;
			new_part1 -= c;
		}

		if (new_part1 < 0 || new_part2 < 0)
			return false;

		SizeType nsite = geometry_.numberOfSites();
		if (SizeType(new_part1) > nsite || SizeType(new_part2) > nsite)
			return false;
		if (new_part1 == 0 && new_part2 == 0)
			return false;
		if (SizeType(new_part1 + new_part2) > nsite)
			return false; // no double occupancy
		newParts.first  = SizeType(new_part1);
		newParts.second = SizeType(new_part2);
		return true;
	}

	void calcDiagonalElements(typename PsimagLite::Vector<RealType>::Type& diag,
	                          const BasisBaseType&                         basis) const
	{
		const RealType zero_point_twenty_five = 0.25;
		SizeType       hilbert                = basis.size();
		SizeType       nsite                  = geometry_.numberOfSites();

		// Calculate diagonal elements
		for (SizeType ispace = 0; ispace < hilbert; ispace++) {
			WordType          ket1 = basis(ispace, SPIN_UP);
			WordType          ket2 = basis(ispace, SPIN_DOWN);
			ComplexOrRealType s    = 0;
			for (SizeType i = 0; i < nsite; i++) {

				int nniup   = 0;
				int nnidown = 0;
				for (SizeType orb = 0; orb < mp_.orbitals; ++orb) {
					nniup += basis.isThereAnElectronAt(
					    ket1, ket2, i, SPIN_UP, orb);
					nnidown += basis.isThereAnElectronAt(
					    ket1, ket2, i, SPIN_DOWN, orb);
				}
				int proi = (nniup + nnidown > 0) ? abs(nniup + nnidown - 1) : 0;

				for (SizeType orb = 0; orb < mp_.orbitals; ++orb) {
					int niup = basis.isThereAnElectronAt(
					    ket1, ket2, i, SPIN_UP, orb);
					int nidown = basis.isThereAnElectronAt(
					    ket1, ket2, i, SPIN_DOWN, orb);

					if (i < mp_.potentialV.size()) {
						s += mp_.potentialV[i + orb * nsite] * niup;
						s += mp_.potentialV[i + orb * nsite
						                    + mp_.orbitals * nsite]
						    * nidown;
					}

					for (SizeType j = i + 1; j < nsite; j++) {

						int nnjup   = 0;
						int nnjdown = 0;
						for (SizeType orb2 = 0; orb2 < mp_.orbitals;
						     ++orb2) {
							nnjup += basis.isThereAnElectronAt(
							    ket1, ket2, j, SPIN_UP, orb2);
							nnjdown += basis.isThereAnElectronAt(
							    ket1, ket2, j, SPIN_DOWN, orb2);
						}

						int proj = (nnjup + nnjdown > 0)
						    ? abs(nnjup + nnjdown - 1)
						    : 0;

						for (SizeType orb2 = 0; orb2 < mp_.orbitals;
						     ++orb2) {
							int njup = basis.isThereAnElectronAt(
							    ket1, ket2, j, SPIN_UP, orb2);
							int njdown = basis.isThereAnElectronAt(
							    ket1, ket2, j, SPIN_DOWN, orb2);
							int proij
							    = (mp_.orbitals > 1) ? proi * proj : 1;
							// Sz Sz term:

							s += proij * (niup - nidown)
							    * (njup - njdown)
							    * jzz_(i + j * nsite,
							           orb + orb2 * mp_.orbitals)
							    * zero_point_twenty_five;
							// ni nj term
							s += proij * (niup + nidown)
							    * (njup + njdown)
							    * w_(i + j * nsite,
							         orb + orb2 * mp_.orbitals);
						}
					}
				}
			}

			assert(fabs(PsimagLite::imag(s)) < 1e-12);
			diag[ispace] = PsimagLite::real(s);
		}
	}

	void setHoppingTerm(SparseRowType&       sparseRow,
	                    const WordType&      ket1,
	                    const WordType&      ket2,
	                    SizeType             i,
	                    SizeType             orb,
	                    const BasisBaseType& basis) const
	{
		WordType s1i = (ket1 & BasisType::bitmask(i * mp_.orbitals + orb));
		if (s1i > 0)
			s1i = 1;
		WordType s2i = (ket2 & BasisType::bitmask(i * mp_.orbitals + orb));
		if (s2i > 0)
			s2i = 1;

		SizeType nsite = geometry_.numberOfSites();

		// Hopping term
		for (SizeType j = 0; j < nsite; j++) {
			for (SizeType orb2 = 0; orb2 < mp_.orbitals; ++orb2) {
				if (j < i)
					continue;
				ComplexOrRealType h
				    = hoppings_(i + j * nsite, orb + orb2 * mp_.orbitals);
				if (PsimagLite::real(h) == 0 && PsimagLite::imag(h) == 0)
					continue;
				WordType s1j = (ket1 & BasisType::bitmask(j * mp_.orbitals + orb2));
				if (s1j > 0)
					s1j = 1;
				WordType s2j = (ket2 & BasisType::bitmask(j * mp_.orbitals + orb2));
				if (s2j > 0)
					s2j = 1;

				if (s1i + s1j == 1 && !(s1j == 0 && s2j > 0)
				    && !(s1j > 0 && s2i > 0)) {
					WordType bra1 = ket1
					    ^ (BasisType::bitmask(i * mp_.orbitals + orb)
					       | BasisType::bitmask(j * mp_.orbitals + orb2));
					SizeType temp = basis.perfectIndex(bra1, ket2);
					RealType extra_sign
					    = (s1i == 1) ? LanczosGlobals::FERMION_SIGN : 1;
					RealType tmp2
					    = basis_.doSign(ket1, ket2, i, orb, j, orb2, SPIN_UP);
					ComplexOrRealType c_temp = h * extra_sign * tmp2;
					sparseRow.add(temp, c_temp);
				}

				if (s2i + s2j == 1 && !(s2j == 0 && s1j > 0)
				    && !(s2j > 0 && s1i > 0)) {
					WordType bra2 = ket2
					    ^ (BasisType::bitmask(i * mp_.orbitals + orb)
					       | BasisType::bitmask(j * mp_.orbitals + orb2));
					SizeType temp = basis.perfectIndex(ket1, bra2);
					RealType extra_sign
					    = (s2i == 1) ? LanczosGlobals::FERMION_SIGN : 1;
					RealType tmp2
					    = basis_.doSign(ket1, ket2, i, orb, j, orb2, SPIN_DOWN);
					ComplexOrRealType c_temp = h * extra_sign * tmp2;
					sparseRow.add(temp, c_temp);
				}
			}
		}
	}

	void setSplusSminus(SparseRowType&       sparseRow,
	                    const WordType&      ket1,
	                    const WordType&      ket2,
	                    SizeType             i,
	                    const BasisBaseType& basis) const
	{
		const RealType zero_point_five = 0.5;
		int            nniup           = 0;
		int            nnidown         = 0;
		for (SizeType orb = 0; orb < mp_.orbitals; ++orb) {
			nniup += basis.isThereAnElectronAt(ket1, ket2, i, SPIN_UP, orb);
			nnidown += basis.isThereAnElectronAt(ket1, ket2, i, SPIN_DOWN, orb);
		}

		int proi = (nniup + nnidown > 0) ? abs(nniup + nnidown - 1) : 0;

		for (SizeType orb = 0; orb < mp_.orbitals; ++orb) {

			WordType s1i = (ket1 & BasisType::bitmask(i * mp_.orbitals + orb));
			if (s1i > 0)
				s1i = 1;
			WordType s2i = (ket2 & BasisType::bitmask(i * mp_.orbitals + orb));
			if (s2i > 0)
				s2i = 1;

			SizeType nsite = geometry_.numberOfSites();

			// Splus Sminus
			for (SizeType j = 0; j < nsite; j++) {
				if (j < i)
					continue;
				int nnjup   = 0;
				int nnjdown = 0;
				for (SizeType orb2 = 0; orb2 < mp_.orbitals; ++orb2) {
					nnjup += basis.isThereAnElectronAt(
					    ket1, ket2, j, SPIN_UP, orb2);
					nnjdown += basis.isThereAnElectronAt(
					    ket1, ket2, j, SPIN_DOWN, orb2);
				}

				int proj = (nnjup + nnjdown > 0) ? abs(nnjup + nnjdown - 1) : 0;

				for (SizeType orb2 = 0; orb2 < mp_.orbitals; ++orb2) {
					ComplexOrRealType h
					    = jpm_(i + j * nsite, orb + orb2 * mp_.orbitals)
					    * zero_point_five;
					if (PsimagLite::real(h) == 0 && PsimagLite::imag(h) == 0)
						continue;
					WordType s1j
					    = (ket1 & BasisType::bitmask(j * mp_.orbitals + orb2));
					if (s1j > 0)
						s1j = 1;
					WordType s2j
					    = (ket2 & BasisType::bitmask(j * mp_.orbitals + orb2));
					if (s2j > 0)
						s2j = 1;

					if (s1i == 1 && s1j == 0 && s2i == 0 && s2j == 1) {
						WordType bra1 = ket1
						    ^ BasisType::bitmask(i * mp_.orbitals + orb);
						bra1 |= BasisType::bitmask(j * mp_.orbitals + orb2);
						WordType bra2 = ket2
						    | BasisType::bitmask(i * mp_.orbitals + orb);
						bra2 ^= BasisType::bitmask(j * mp_.orbitals + orb2);
						SizeType temp = basis.perfectIndex(bra1, bra2);
						int proij = (mp_.orbitals > 1) ? proi * proj : 1;
						sparseRow.add(
						    temp,
						    proij * h
						        * signSplusSminus(i * mp_.orbitals + orb,
						                          j * mp_.orbitals + orb2,
						                          bra1,
						                          bra2));
					}

					if (s1i == 0 && s1j == 1 && s2i == 1 && s2j == 0) {
						WordType bra1 = ket1
						    | BasisType::bitmask(i * mp_.orbitals + orb);
						bra1 ^= BasisType::bitmask(j * mp_.orbitals + orb2);
						WordType bra2 = ket2
						    ^ BasisType::bitmask(i * mp_.orbitals + orb);
						bra2 |= BasisType::bitmask(j * mp_.orbitals + orb2);
						SizeType temp = basis.perfectIndex(bra1, bra2);
						int proij = (mp_.orbitals > 1) ? proi * proj : 1;
						sparseRow.add(
						    temp,
						    proij * h
						        * signSplusSminus(i * mp_.orbitals + orb,
						                          j * mp_.orbitals + orb2,
						                          bra1,
						                          bra2));
					}
				}
			}
		}
	}

	RealType
	signSplusSminus(SizeType i, SizeType j, const WordType& bra1, const WordType& bra2) const
	{
		int s = 1;
		if (j > 0)
			s *= parityFrom(0, j - 1, bra2);
		if (i > 0)
			s *= parityFrom(0, i - 1, bra2);
		if (i > 0)
			s *= parityFrom(0, i - 1, bra1);
		if (j > 0)
			s *= parityFrom(0, j - 1, bra1);
		return s;
	}

	// from i to j including i and j
	// assumes i<=j
	int parityFrom(SizeType i, SizeType j, const WordType& ket) const
	{
		if (i == j)
			return (BasisType::bitmask(j) & ket) ? -1 : 1;
		assert(i < j);
		// j>i>=0 now
		WordType mask = ket;
		mask &= ((1 << (i + 1)) - 1) ^ ((1 << j) - 1);
		int s = (PsimagLite::BitManip::count(mask) & 1) ? -1 : 1;
		// Is there something of this species at i?
		if (BasisType::bitmask(i) & ket)
			s = -s;
		// Is there something of this species at j?
		if (BasisType::bitmask(j) & ket)
			s = -s;
		return s;
	}

	const ParametersModelType                             mp_;
	const GeometryType&                                   geometry_;
	BasisType                                             basis_;
	PsimagLite::Matrix<ComplexOrRealType>                 hoppings_;
	PsimagLite::Matrix<ComplexOrRealType>                 jpm_;
	PsimagLite::Matrix<ComplexOrRealType>                 jzz_;
	PsimagLite::Matrix<ComplexOrRealType>                 w_;
	mutable typename PsimagLite::Vector<BasisType*>::Type garbage_;
}; // class TjMultiOrb
} // namespace LanczosPlusPlus
#endif
