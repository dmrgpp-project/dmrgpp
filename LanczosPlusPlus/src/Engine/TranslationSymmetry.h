/*
Copyright (c) 2009-2014, UT-Battelle, LLC
All rights reserved

[Lanczos++, Version 1.0.0]

*********************************************************
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED.

Please see full open source license included in file LICENSE.
*********************************************************

*/

#ifndef TRANSLATION_SYMM_H
#define TRANSLATION_SYMM_H
#include "LanczosGlobals.h"
#include <PsimagLite/CrsMatrix.h>
#include <PsimagLite/ProgressIndicator.h>
#include <PsimagLite/SparseVector.h>
#include <PsimagLite/Vector.h>
#include <iostream>

namespace LanczosPlusPlus {

struct TranslationItem {
	SizeType type;
	SizeType i;
};

template <typename RealType> class Kspace {

public:

	Kspace(SizeType len)
	    : data_(len, len)
	    , blockSizes_(len, 0)
	{

		for (SizeType i = 0; i < len - 1; i++) {
			data_(i, i + 1) = 1;
		}
		data_(len - 1, 0) = 1;
		typename PsimagLite::Vector<RealType>::Type e(len);
		diag(data_, e, 'V');
	}

	SizeType size() const { return data_.n_row(); }

	const RealType& operator()(SizeType i, SizeType j) const { return data_(i, j); }

	void setBlockSize(SizeType k, SizeType s)
	{
		std::cout << "BLOCKSZIZE[" << k << "]=" << s << "\n";
		blockSizes_[k] = s;
	}

	SizeType blockSizes(SizeType k) const { return blockSizes_[k]; }

	SizeType blockSize() const
	{
		SizeType sum = 0;
		for (SizeType i = 0; i < blockSizes_.size(); i++)
			sum += blockSizes_[i];
		return sum;
	}

private:

	PsimagLite::Matrix<RealType>       data_;
	PsimagLite::Vector<SizeType>::Type blockSizes_;
};

template <typename GeometryType, typename BasisType, typename KspaceType>
class ClassRepresentatives {

	typedef typename BasisType::WordType WordType;
	typedef TranslationItem              ItemType;

public:

	ClassRepresentatives(const BasisType&    basis,
	                     const GeometryType& geometry,
	                     const KspaceType&   kspace)
	    : basis_(basis)
	    , geometry_(geometry)
	    , data_(basis.size())
	{
		SizeType                       hilbert = basis.size();
		PsimagLite::Vector<bool>::Type seen(hilbert, false);

		for (SizeType ispace = 0; ispace < hilbert; ispace++) {

			for (SizeType k = 0; k < kspace.size(); k++) {
				typename PsimagLite::Vector<WordType>::Type y
				    = translateInternal(ispace, k);
				SizeType y_index = basis.perfectIndex(y);
				if (!seen[y_index]) {
					seen[y_index]      = true;
					data_[yIndex].type = k;
					data_[yIndex].i    = ispace;
				}
			}
		}
	}

	SizeType translate(SizeType state, SizeType k) const
	{
		for (SizeType i = 0; i < data_.size(); i++) {
			if (data_[i].i == state && data_[i].type == k)
				return i;
		}

		return state;
	}

private:

	typename PsimagLite::Vector<WordType>::Type translateInternal(SizeType state,
	                                                              SizeType k) const
	{
		SizeType                                    number_of_dofs = basis_.dofs();
		typename PsimagLite::Vector<WordType>::Type y(number_of_dofs);

		for (SizeType dof = 0; dof < number_of_dofs; dof++) {
			WordType x = basis_(state, dof);
			y[dof]     = translateInternal2(x, k);
		}
		std::cout << "translation of " << basis_(state, 0) << " " << basis_(state, 1);
		std::cout << "is " << y[0] << " " << y[1] << "\n";
		return y;
	}

	WordType translateInternal2(WordType state, SizeType k) const
	{
		SizeType number_of_sites = geometry_.numberOfSites();
		SizeType term_id         = 0;
		WordType x               = state;
		WordType y               = 0;
		SizeType diry            = 1;
		for (SizeType site = 0; site < number_of_sites; site++) {
			SizeType t_site = geometry_.term(term_id).translate(site, diry, k);
			SizeType this_site_content = x & 1;
			x >>= 1; // go to next site
			addTo(y, this_site_content, t_site);
			if (!x)
				break;
		}

		return y;
	}

	void addTo(WordType& yy, SizeType what, SizeType site) const
	{
		if (what == 0)
			return;
		WordType mask = (1 << site);
		yy |= mask;
	}

	const BasisType&                          basis_;
	const GeometryType&                       geometry_;
	PsimagLite::Vector<TranslationItem>::Type data_;
};

template <typename GeometryType_, typename BasisType> class TranslationSymmetry {

	typedef typename GeometryType_::ComplexOrRealType                  ComplexOrRealType;
	typedef PsimagLite::Matrix<ComplexOrRealType>                      MatrixType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type         RealType;
	typedef LanczosGlobals::WordType                                   WordType;
	typedef PsimagLite::SparseVector<ComplexOrRealType>                SparseVectorType;
	typedef Kspace<RealType>                                           KspaceType;
	typedef ClassRepresentatives<GeometryType_, BasisType, KspaceType> ClassRepresentativesType;
	typedef std::pair<PsimagLite::Vector<SizeType>::Type, SizeType>    BufferItemType;
	typedef typename PsimagLite::Vector<SparseVectorType>::Type        VectorSparseVectorType;

public:

	typedef GeometryType_                                        GeometryType;
	typedef PsimagLite::CrsMatrix<ComplexOrRealType>             SparseMatrixType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
	typedef typename PsimagLite::Vector<RealType>::Type          VectorRealType;
	typedef typename PsimagLite::Vector<VectorType>::Type        VectorVectorType;

	TranslationSymmetry(const BasisType&    basis,
	                    const GeometryType& geometry,
	                    PsimagLite::String  options)
	    : progress_("TranslationSymmetry")
	    , transform_(basis.size(), basis.size())
	    , kspace_(geometry.length(1, 0))
	    , matrixStored_(kspace_.size())
	    , pointer_(0)
	    , printMatrix_(options.find("printmatrix") != PsimagLite::String::npos)
	{
		ClassRepresentativesType reps(basis, geometry, kspace_);

		SizeType               hilbert = basis.size();
		VectorSparseVectorType bag;
		for (SizeType k = 0; k < kspace_.size(); k++) {
			SizeType block_size = 0;
			for (SizeType ispace = 0; ispace < hilbert; ispace++) {
				typename PsimagLite::Vector<ComplexOrRealType>::Type v(hilbert);
				eikrTr(v, ispace, k, reps);
				SparseVectorType sparse_v(v);
				sparse_v.sort();
				if (!checkForOrthogonality(sparse_v, bag))
					continue;
				bag.push_back(sparse_v);
				block_size++;
			}
			kspace_.setBlockSize(k, block_size);
		}

		if (kspace_.blockSize() != hilbert) {
			std::cout << "Blocksizes summed=" << kspace_.blockSize();
			std::cout << " but hilbert=" << hilbert << "\n";
			throw std::runtime_error("error!\n");
		}
		setTransform(bag);
		//			checkTransform();
	}

	template <typename SomeModelType>
	void init(const SomeModelType& model, const BasisType& basis)
	{
		SparseMatrixType matrix2;
		model.setupHamiltonian(matrix2, basis);
		transformMatrix(matrixStored_, matrix2);

		if (matrixStored_.size() == 0)
			return;
		int nrows = matrixStored_[0].rows();
		if (printMatrix_) {
			if (nrows > 40)
				throw PsimagLite::RuntimeError("printMatrix too big\n");
			std::cout << matrixStored_[0].toDense();
		}
	}

	SizeType rows() const { return matrixStored_[pointer_].rows(); }

	template <typename SomeVectorType>
	void matrixVectorProduct(SomeVectorType& x, SomeVectorType const& y) const
	{
		return matrixStored_[pointer_].matrixVectorProduct(x, y);
	}

	void transformMatrix(typename PsimagLite::Vector<SparseMatrixType>::Type& matrix1,
	                     const SparseMatrixType&                              matrix) const
	{
		SparseMatrixType r_t;
		transposeConjugate(r_t, transform_);

		if (matrix.rows() < 40)
			printFullMatrix(matrix, "originalHam");
		SparseMatrixType tmp;
		multiply(tmp, matrix, r_t);

		SparseMatrixType matrix2;
		multiply(matrix2, transform_, tmp);

		if (matrix2.rows() < 40)
			printFullMatrix(matrix2, "HamiltonianTransformed");
		matrix1.clear();
		split(matrix1, matrix2);
	}

	void transform(VectorVectorType& zs, SizeType offset)
	{
		VectorType gstmp(transform_.rows(), 0);

		const SizeType excited_plus_one = zs.size();

		for (SizeType i = 0; i < excited_plus_one; ++i) {
			LanczosGlobals::transform(zs[i], offset, gstmp, transform_);
		}
	}

	SizeType sectors() const { return kspace_.size(); }

	void setPointer(SizeType p) { pointer_ = p; }

	PsimagLite::String name() const { return "translation"; }

	void fullDiag(VectorRealType& eigs, MatrixType& fm) const
	{
		if (matrixStored_[pointer_].rows() > 1000)
			throw PsimagLite::RuntimeError("fullDiag too big\n");

		fm = matrixStored_[pointer_].toDense();
		diag(fm, eigs, 'V');

		if (!printMatrix_)
			return;

		for (SizeType i = 0; i < eigs.size(); i++)
			std::cout << eigs[i] << "\n";
		std::cout << fm;
	}

private:

	void addTo(WordType& yy, SizeType what, SizeType site) const
	{
		if (what == 0)
			return;
		WordType mask = (1 << site);
		yy |= mask;
	}

	void setTransform(const VectorSparseVectorType& bag)
	{
		SizeType counter = 0;
		for (SizeType row = 0; row < bag.size(); row++) {
			transform_.setRow(row, counter);
			for (SizeType k = 0; k < bag[row].indices(); k++) {
				transform_.pushCol(bag[row].index(k));
				transform_.pushValue(bag[row].value(k));
				counter++;
			}
		}

		transform_.setRow(transform_.rows(), counter);
		transform_.checkValidity();
		if (transform_.rows() < 40)
			printFullMatrix(transform_, "transform");
	}

	void eikrTr(typename PsimagLite::Vector<std::complex<RealType>>::Type& v,
	            SizeType                                                   ispace,
	            SizeType                                                   k,
	            const ClassRepresentativesType&                            reps) const
	{
		for (SizeType r = 0; r < kspace_.size(); r++) {
			SizeType jspace = reps.translate(ispace, r);
			RealType tmp    = 2 * M_PI * k * r / RealType(kspace_.size());
			v[jspace]       = ComplexOrRealType(cos(tmp), sin(tmp));
		}
	}

	void eikrTr(typename PsimagLite::Vector<RealType>::Type&,
	            SizeType,
	            SizeType,
	            const ClassRepresentativesType&) const
	{
		throw PsimagLite::RuntimeError("eikrTr: not for real template\n");
	}

	bool checkForOrthogonality(const SparseVectorType&       sparseV,
	                           const VectorSparseVectorType& bag) const
	{
		for (SizeType i = 0; i < bag.size(); i++) {
			ComplexOrRealType sp = sparseV.scalarProduct(bag[i]);
			if (PsimagLite::norm(sp) > 1e-8)
				return false;
		}
		return true;
	}

	void split(typename PsimagLite::Vector<SparseMatrixType>::Type& matrix,
	           const SparseMatrixType&                              matrix2) const
	{
		SizeType offset = 0;
		for (SizeType i = 0; i < kspace_.size(); i++) {
			SizeType block_size = kspace_.blockSizes(i);
			if (block_size == 0)
				continue;
			std::cout << "BLOCKSIZE=" << block_size << "\n";
			SparseMatrixType m(block_size, block_size);
			SizeType         counter = 0;
			for (SizeType row = 0; row < block_size; row++) {
				m.setRow(row, counter);
				SizeType global_row = row + offset;
				SizeType start      = matrix2.getRowPtr(global_row);
				SizeType end        = matrix2.getRowPtr(global_row + 1);
				for (SizeType k = start; k < end; k++) {
					ComplexOrRealType val = matrix2.getValue(k);
					if (PsimagLite::norm(val) < 1e-8)
						continue;
					SizeType global_col = matrix2.getCol(k);

					assert(globalCol >= offset);
					SizeType col = global_col - offset;

					assert(col < blockSize);
					m.pushCol(col);
					m.pushValue(val);
					counter++;
				}
			}
			m.setRow(block_size, counter);
			m.checkValidity();
			matrix.push_back(m);
			offset += block_size;
		}
	}

	void transform(VectorType& gs, SizeType offset, VectorType& gstmp)
	{
		for (SizeType i = 0; i < gs.size(); i++) {
			assert(i + offset < gstmp.size());
			gstmp[i + offset] = gs[i];
		}

		SparseMatrixType r_t;
		transposeConjugate(r_t, transform_);
		gs.clear();
		gs.resize(transform_.rows());
		multiply(gs, r_t, gstmp);
	}

	PsimagLite::ProgressIndicator                       progress_;
	SparseMatrixType                                    transform_;
	KspaceType                                          kspace_;
	typename PsimagLite::Vector<SparseMatrixType>::Type matrixStored_;
	SizeType                                            pointer_;
	bool                                                printMatrix_;
}; // class TranslationSymmetry
} // namespace Dmrg

#endif // TRANSLATION_SYMM_H
