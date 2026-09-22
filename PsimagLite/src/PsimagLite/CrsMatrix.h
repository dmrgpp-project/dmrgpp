/*
Copyright (c) 2009-2012-2018, UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 2.]
[by G.A., Oak Ridge National Laboratory]

UT Battelle Open Source Software License 11242008

OPEN SOURCE LICENSE

Subject to the conditions of this License, each
contributor to this software hereby grants, free of
charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), a
perpetual, worldwide, non-exclusive, no-charge,
royalty-free, irrevocable copyright license to use, copy,
modify, merge, publish, distribute, and/or sublicense
copies of the Software.

1. Redistributions of Software must retain the above
copyright and license notices, this list of conditions,
and the following disclaimer.  Changes or modifications
to, or derivative works of, the Software should be noted
with comments and the contributor and organization's
name.

2. Neither the names of UT-Battelle, LLC or the
Department of Energy nor the names of the Software
contributors may be used to endorse or promote products
derived from this software without specific prior written
permission of UT-Battelle.

3. The software and the end-user documentation included
with the redistribution, with or without modification,
must include the following acknowledgment:

"This product includes software produced by UT-Battelle,
LLC under Contract No. DE-AC05-00OR22725  with the
Department of Energy."

*********************************************************
DISCLAIMER

THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT OWNER, CONTRIBUTORS, UNITED STATES GOVERNMENT,
OR THE UNITED STATES DEPARTMENT OF ENERGY BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED
STATES DEPARTMENT OF ENERGY, NOR THE COPYRIGHT OWNER, NOR
ANY OF THEIR EMPLOYEES, REPRESENTS THAT THE USE OF ANY
INFORMATION, DATA, APPARATUS, PRODUCT, OR PROCESS
DISCLOSED WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*********************************************************

*/
/** \ingroup DMRG */
/*@{*/

/*! \file CrsMatrix.h
 *
 *  A class to represent a sparse matrix in Compressed Row Storage
 *
 */

#ifndef CRSMATRIX_HEADER_H
#define CRSMATRIX_HEADER_H
#include "BLAS.h"
#include "Complex.h"
#include "Io/IoSerializerStub.h"
#include "Matrix.h"
#include "NotMpi.h"
#include "Sort.h"
#include "loki/TypeTraits.h"
#include <algorithm>
#include <cassert>
#include <fstream>
#include <unordered_map>
#include <vector>

namespace PsimagLite {

//! A Sparse Matrix in Compressed Row Storage (CRS) format.
/**
        The CRS format puts the subsequent nonzero elements of the matrix rows
        in contiguous memory locations. We create 3 vectors: one for complex
   numbers containing the values of the matrix entries and the other two for
   integers ($colind$ and $rowptr$). The vector $values$ stores the values of
   the non-zero elements of the matrix, as they are traversed in a row-wise
   fashion. The $colind$ vector stores the column indices of the elements of the
   $values$ vector. That is, if $values[k] = a[i][j]$ then $colind[k] = j$. The
   $rowptr$ vector stores the locations in the $values$ vector that start a row,
   that is $values[k] = a[i][j]$ if $rowptr[i] \le i < rowptr[i + 1]$. By
   convention, we define $rowptr[N_{dim}]$ to be equal to the number of non-zero
   elements, $n_z$, in the matrix. The storage savings of this approach are
   significant because instead of
        storing $N_{dim}^2$ elements, we need only $2n_z + N_{dim} + 1$ storage
   locations.\\ To illustrate how the CRS format works, consider the
   non-symmetric matrix defined by \begin{equation}
                A=\left[\begin{tabular}{llllll}

                10 &  0 & 0 & 0  & -2 & 0 \\
                3 &  9 &  0 &  0 &  0 &  3 \\
                0 &  7 &  8 &  7 &  0 &  0 \\
                3 &  0 &  8 &  7  & 5 &  0 \\
                0 &   8 &  0 &  9 &  9 & 13 \\
                0 &  4 &  0 &  0 &  2&  -1 \\
        \end{tabular}\right]\end{equation}
        The CRS format for this matrix is then specified by the arrays:\\
        \begin{tt}
                values = [10 -2  3  9  3  7  8  7  3 ... 9 13  4  2 -1 ]\\
                colind = [ 0  4  0  1  5  1  2  3  0 ... 4  5  1  4  5 ]\\
                rowptr = [ 0  2  5  8 12 16 19 ]\\
        \end{tt}
        */
template <class T> class CrsMatrix {

public:

	using MatrixElementType = T;
	using ValueType         = T;

	CrsMatrix()
	    : nrow_(0)
	    , ncol_(0)
	{ }

	~CrsMatrix() { }

	CrsMatrix(SizeType nrow, SizeType ncol)
	    : nrow_(nrow)
	    , ncol_(ncol)
	{
		resize(nrow, ncol);
	}

	CrsMatrix(SizeType nrow, SizeType ncol, SizeType nonzero)
	    : nrow_(nrow)
	    , ncol_(ncol)
	{
		resize(nrow, ncol, nonzero);
	}

	template <typename S> CrsMatrix(const CrsMatrix<S>& a)
	{
		colind_ = a.colind_;
		rowptr_ = a.rowptr_;
		values_ = a.values_;
		nrow_   = a.nrow_;
		ncol_   = a.ncol_;
	}

	template <typename S> CrsMatrix(const CrsMatrix<std::complex<S>>& a)
	{
		colind_ = a.colind_;
		rowptr_ = a.rowptr_;
		values_ = a.values_;
		nrow_   = a.nrow_;
		ncol_   = a.ncol_;
	}

	explicit CrsMatrix(const Matrix<T>& a)
	{
		int    counter = 0;
		double eps     = 0;

		resize(a.rows(), a.cols());

		for (SizeType i = 0; i < a.rows(); i++) {
			setRow(i, counter);
			for (SizeType j = 0; j < a.cols(); j++) {
				if (PsimagLite::norm(a(i, j)) <= eps)
					continue;
				pushValue(a(i, j));
				pushCol(j);
				++counter;
			}
		}

		setRow(a.rows(), counter);
	}

	CrsMatrix(SizeType                        rank, // square matrix ONLY for now
	          const Vector<SizeType>::Type&   rows2,
	          const Vector<SizeType>::Type&   cols,
	          const typename Vector<T>::Type& vals)
	    : rowptr_(rank + 1)
	    , nrow_(rank)
	    , ncol_(rank)
	{
		Sort<Vector<SizeType>::Type> s;
		Vector<SizeType>::Type       iperm(rows2.size());
		Vector<SizeType>::Type       rows = rows2;
		s.sort(rows, iperm);
		SizeType counter  = 0;
		SizeType prev_row = rows[0] + 1;
		for (SizeType i = 0; i < rows.size(); i++) {
			SizeType row = rows[i];
			if (prev_row != row) {
				// add new row
				rowptr_[row] = counter++;
				prev_row     = row;
			}

			colind_.push_back(cols[iperm[i]]);
			values_.push_back(vals[iperm[i]]);
		}

		SizeType last_non_zero_row = rows[rows.size() - 1];
		for (SizeType i = last_non_zero_row + 1; i <= rank; ++i)
			rowptr_[i] = counter;
	}

	// start closure ctors

	CrsMatrix(
	    const std::ClosureOperator<CrsMatrix, CrsMatrix, std::ClosureOperations::OP_MULT>& c)
	{
		CrsMatrix&       x = *this;
		const CrsMatrix& y = c.r1;
		const CrsMatrix& z = c.r2;
		multiply(x, y, z);
	}

	CrsMatrix(const std::ClosureOperator<T, CrsMatrix, std::ClosureOperations::OP_MULT>& c)
	{
		*this = c.r2;
		this->values_ *= c.r1;
	}

	// end all ctors

	const CrsMatrix<T>& toCRS() const { return *this; }

	void resize(SizeType nrow, SizeType ncol)
	{
		colind_.clear();
		values_.clear();
		rowptr_.clear();
		rowptr_.resize(nrow + 1);
		nrow_ = nrow;
		ncol_ = ncol;
	}

	void clear()
	{
		colind_.clear();
		values_.clear();
		rowptr_.clear();
		nrow_ = ncol_ = 0;
	}

	void resize(SizeType nrow, SizeType ncol, SizeType nonzero)
	{
		nrow_ = nrow;
		ncol_ = ncol;

		// ------------------------------------
		// Note arrays are not cleared out
		// arrays should retain original values
		// ------------------------------------
		rowptr_.resize(nrow_ + 1);
		colind_.resize(nonzero);
		values_.resize(nonzero);
	}

	void reserve(SizeType nonzero)
	{
		// -------------------------------------------------
		// increase internal capacity
		// to avoid repeated allocation expansion and copies
		// -------------------------------------------------
		colind_.reserve(nonzero);
		values_.reserve(nonzero);
	}

	void setRow(SizeType n, SizeType v)
	{
		assert(n < rowptr_.size());
		rowptr_[n] = v;
	}

	void setCol(int n, int v) { colind_[n] = v; }

	void setColCheck(int n, int v)
	{
		if (((size_t)n) == (colind_.size() + 1)) {
			colind_.push_back(v);
		} else {
			colind_[n] = v;
		};
	}

	void setValues(int n, const T& v) { values_[n] = v; }

	void setValuesCheck(int n, const T& v)
	{
		if (((size_t)n) == (values_.size() + 1)) {
			values_.push_back(v);
		} else {
			values_[n] = v;
		};
	}

	void operator*=(T x) { values_ *= x; }

	bool operator==(const CrsMatrix<T>& op) const
	{
		return (nrow_ == op.nrow_ && ncol_ == op.ncol_ && rowptr_ == op.rowptr_
		        && colind_ == op.colind_ && values_ == op.values_);
	}

	template <typename VerySparseMatrixType>
	typename EnableIf<!std::IsClosureLike<VerySparseMatrixType>::True, void>::Type
	operator=(const VerySparseMatrixType& m)
	{
		if (!m.sorted())
			throw RuntimeError("CrsMatrix: VerySparseMatrix must be sorted\n");

		clear();
		SizeType non_zeros = m.nonZeros();
		resize(m.rows(), m.cols(), non_zeros);

		SizeType counter = 0;
		for (SizeType i = 0; i < m.rows(); ++i) {
			setRow(i, counter);

			while (counter < non_zeros && m.getRow(counter) == i) {
				colind_[counter] = m.getColumn(counter);
				values_[counter] = m.getValue(counter);
				counter++;
			}
		}

		setRow(m.rows(), counter);
		checkValidity();
	}

	void operator+=(const CrsMatrix& m)
	{
		CrsMatrix                           c;
		static const typename Real<T>::Type f1 = 1.0;
		add(c, m, f1);
		*this = c;
	}

	SizeType nonZeros() const
	{
		if (nrow_ >= 1) {
			assert(rowptr_.size() == 1 + nrow_);

			assert(static_cast<SizeType>(rowptr_[nrow_]) == colind_.size());
			assert(static_cast<SizeType>(rowptr_[nrow_]) == values_.size());

			return colind_.size();
		} else {
			return 0;
		};
	}

	/** performs x = x + A * y
	 ** where x and y are vectors and A is a sparse matrix in
	 ** row-compressed format */
	template <typename VectorLikeType>
	void matrixVectorProduct(VectorLikeType& x, const VectorLikeType& y) const
	{
		assert(x.size() == y.size());
		for (SizeType i = 0; i < y.size(); i++) {
			assert(i + 1 < rowptr_.size());
			for (int j = rowptr_[i]; j < rowptr_[i + 1]; j++) {
				assert(SizeType(j) < values_.size());
				assert(SizeType(j) < colind_.size());
				assert(SizeType(colind_[j]) < y.size());
				x[i] += values_[j] * y[colind_[j]];
			}
		}
	}

#ifndef NO_DEPRECATED_ALLOWED
	int nonZero() const { return colind_.size(); } // DEPRECATED, use nonZeros()
#endif

	SizeType rows() const { return nrow_; }

	SizeType cols() const { return ncol_; }

	void swap(CrsMatrix& other)
	{
		this->rowptr_.swap(other.rowptr_);
		this->colind_.swap(other.colind_);
		this->values_.swap(other.values_);
		SizeType nrow = this->nrow_;
		this->nrow_   = other.nrow_;
		other.nrow_   = nrow;
		SizeType ncol = this->ncol_;
		this->ncol_   = other.ncol_;
		other.ncol_   = ncol;
	}

	void pushCol(SizeType i) { colind_.push_back(i); }

	void pushValue(T const& value) { values_.push_back(value); }

	//! Make a diagonal CRS matrix with value "value"
	void makeDiagonal(SizeType row, T const& value = 0)
	{
		nrow_ = row;
		ncol_ = row;
		rowptr_.resize(row + 1);
		values_.resize(row);
		colind_.resize(row);

		for (SizeType i = 0; i < row; i++) {
			values_[i] = value;
			colind_[i] = i;
			rowptr_[i] = i;
		}

		rowptr_[row] = row;
	}

	const int& getRowPtr(SizeType i) const
	{
		assert(i < rowptr_.size());
		return rowptr_[i];
	}

	const int& getCol(SizeType i) const
	{
		assert(i < colind_.size());
		return colind_[i];
	}

	const T& getValue(SizeType i) const
	{
		assert(i < values_.size());
		return values_[i];
	}

	void conjugate()
	{
		SizeType n = values_.size();
		for (SizeType i = 0; i < n; ++i)
			values_[i] = PsimagLite::conj(values_[i]);
	}

	Matrix<T> toDense() const
	{
		Matrix<T> m;
		crsMatrixToFullMatrix(m, *this);
		return m;
	}

	void checkValidity() const
	{
#ifndef NDEBUG
		SizeType n = nrow_;
		assert(n + 1 == rowptr_.size());
		assert(static_cast<SizeType>(rowptr_[nrow_]) == colind_.size());
		assert(values_.size() == colind_.size());
		assert(nrow_ > 0 && ncol_ > 0);
		typename Vector<SizeType>::Type p(ncol_, 0);
		for (SizeType i = 0; i < n; i++) {
			assert(rowptr_[i] <= rowptr_[i + 1]);
			for (int k = rowptr_[i]; k < rowptr_[i + 1]; k++) {
				SizeType col = colind_[k];
				assert(col < p.size());
				assert(p[col] == 0);
				p[col] = 1;
			}

			for (int k = rowptr_[i]; k < rowptr_[i + 1]; k++)
				p[colind_[k]] = 0;
		}
#endif
	}

	// closures operators start

	template <typename T1>
	typename EnableIf<Loki::TypeTraits<T1>::isArith || IsComplexNumber<T1>::True,
	                  CrsMatrix>::Type
	operator=(const std::ClosureOperator<T1, CrsMatrix, std::ClosureOperations::OP_MULT>& c)
	{
		*this = c.r2;
		this->values_ *= c.r1;
		return *this;
	}

	template <typename T1>
	typename EnableIf<Loki::TypeTraits<T1>::isArith || IsComplexNumber<T1>::True,
	                  CrsMatrix>::Type
	operator+=(const std::ClosureOperator<T1, CrsMatrix, std::ClosureOperations::OP_MULT>& c)
	{
		CrsMatrix s;
		add(s, c.r2, c.r1);
		this->swap(s);
		return *this;
	}

	template <typename T1>
	typename EnableIf<Loki::TypeTraits<T1>::isArith || IsComplexNumber<T1>::True,
	                  CrsMatrix>::Type
	operator+=(const std::ClosureOperator<
	           std::ClosureOperator<T1, CrsMatrix, std::ClosureOperations::OP_MULT>,
	           CrsMatrix,
	           std::ClosureOperations::OP_MULT>& c)
	{
		CrsMatrix s;
		multiply(s, c.r1.r2, c.r2);
		CrsMatrix s2;
		add(s2, s, c.r1.r1);
		this->swap(s2);
		return *this;
	}

	// closures operators end

	void send(int root, int tag, MPI::CommType mpiComm)
	{
		MPI::send(nrow_, root, tag, mpiComm);
		MPI::send(ncol_, root, tag + 1, mpiComm);
		MPI::send(rowptr_, root, tag + 2, mpiComm);
		MPI::send(colind_, root, tag + 3, mpiComm);
		MPI::send(values_, root, tag + 4, mpiComm);
	}

	void recv(int root, int tag, MPI::CommType mpiComm)
	{
		MPI::recv(nrow_, root, tag, mpiComm);
		MPI::recv(ncol_, root, tag + 1, mpiComm);
		MPI::recv(rowptr_, root, tag + 2, mpiComm);
		MPI::recv(colind_, root, tag + 3, mpiComm);
		MPI::recv(values_, root, tag + 4, mpiComm);
	}

	void write(String                  label,
	           IoSerializer&           ioSerializer,
	           IoSerializer::WriteMode mode = IoSerializer::NO_OVERWRITE) const
	{
		if (nrow_ > 0)
			checkValidity();
		if (mode != IoSerializer::ALLOW_OVERWRITE)
			ioSerializer.createGroup(label);
		ioSerializer.write(label + "/nrow_", nrow_, mode);
		ioSerializer.write(label + "/ncol_", ncol_, mode);
		if (nrow_ == 0 || ncol_ == 0)
			return;
		ioSerializer.write(label + "/rowptr_", rowptr_, mode);
		assert(rowptr_.size() == nrow_ + 1);
		if (rowptr_[nrow_] == 0)
			return;
		ioSerializer.write(label + "/colind_", colind_, mode);
		ioSerializer.write(label + "/values_", values_, mode);
	}

	void overwrite(String label, IoSerializer& ioSerializer) const
	{
		write(label, ioSerializer, IoSerializer::ALLOW_OVERWRITE);
	}
	void read(String label, IoSerializer& ioSerializer)
	{
		ioSerializer.read(nrow_, label + "/nrow_");
		ioSerializer.read(ncol_, label + "/ncol_");
		if (nrow_ == 0 || ncol_ == 0)
			return;
		ioSerializer.read(rowptr_, label + "/rowptr_");
		assert(rowptr_.size() == nrow_ + 1);
		if (rowptr_[nrow_] == 0)
			return;
		ioSerializer.read(colind_, label + "/colind_");
		ioSerializer.read(values_, label + "/values_");
		checkValidity();
	}

	template <typename S> friend bool isZero(const CrsMatrix<S>&, double);

	template <typename S> friend typename Real<S>::Type norm2(const CrsMatrix<S>& m);

	template <typename S>
	friend std::ostream& operator<<(std::ostream& os, const CrsMatrix<S>& m);

	template <class S> friend void difference(const CrsMatrix<S>& A, const CrsMatrix<S>& B);

	template <typename S> friend void mpiBroadcast(CrsMatrix<S>* v, int rank);

	template <typename S> friend void mpiSend(CrsMatrix<S>* v, int iproc, int i);

	template <typename S> friend void mpiRecv(CrsMatrix<S>* v, int iproc, int i);

	template <typename CrsMatrixType>
	friend std::istream& operator>>(std::istream& is, CrsMatrix<CrsMatrixType>& m);

	template <typename S> friend void bcast(CrsMatrix<S>& m);

private:

	template <typename T1>
	typename std::enable_if<std::is_same<T1, typename Real<T>::Type>::value
	                            || std::is_same<T1, T>::value,
	                        void>::type
	add(CrsMatrix<T>& c, const CrsMatrix<T>& m, const T1& t1) const
	{
		assert(m.rows() == m.cols());
		const T1 one = 1.0;
		if (nrow_ >= m.rows())
			operatorPlus(c, *this, one, m, t1);
		else
			operatorPlus(c, m, t1, *this, one);
	}

	typename Vector<int>::Type rowptr_;
	typename Vector<int>::Type colind_;
	typename Vector<T>::Type   values_;
	SizeType                   nrow_;
	SizeType                   ncol_;
}; // class CrsMatrix

// Companion functions below:

template <typename T> std::ostream& operator<<(std::ostream& os, const CrsMatrix<T>& m)
{
	SizeType n = m.rows();
	if (n == 0) {
		os << "0 0\n";
		return os;
	}

	os << n << " " << m.cols() << "\n";
	for (SizeType i = 0; i < n + 1; i++)
		os << m.rowptr_[i] << " ";
	os << "\n";

	SizeType nonzero = m.nonZeros();
	os << nonzero << "\n";
	for (SizeType i = 0; i < nonzero; i++)
		os << m.colind_[i] << " ";
	os << "\n";

	os << nonzero << "\n";
	for (SizeType i = 0; i < nonzero; i++)
		os << m.values_[i] << " ";
	os << "\n";

	return os;
}

template <typename T> std::istream& operator>>(std::istream& is, CrsMatrix<T>& m)
{
	int n;
	is >> n;
	if (n < 0)
		throw RuntimeError("is>>CrsMatrix(...): Rows must be positive\n");

	int ncol = 0;
	if (ncol < 0)
		throw RuntimeError("is>>CrsMatrix(...): Cols must be positive\n");
	is >> ncol;

	if (n == 0 || ncol == 0)
		return is;

	m.resize(n, ncol);
	for (SizeType i = 0; i < m.rowptr_.size(); i++)
		is >> m.rowptr_[i];

	SizeType nonzero;
	is >> nonzero;
	m.colind_.resize(nonzero);
	for (SizeType i = 0; i < m.colind_.size(); i++)
		is >> m.colind_[i];

	is >> nonzero;
	m.values_.resize(nonzero);
	for (SizeType i = 0; i < m.values_.size(); i++)
		is >> m.values_[i];

	return is;
}

template <typename T> class IsMatrixLike<CrsMatrix<T>> {
public:

	enum
	{
		TRUE = true
	};
};

template <typename S> void bcast(CrsMatrix<S>& m)
{
	MPI::bcast(m.rowptr_);
	MPI::bcast(m.colind_);
	MPI::bcast(m.values_);
	MPI::bcast(m.nrow_);
	MPI::bcast(m.ncol_);
}

//! Transforms a Compressed-Row-Storage (CRS) into a full Matrix (Fast version)
template <typename T> void crsMatrixToFullMatrix(Matrix<T>& m, const CrsMatrix<T>& crsMatrix)
{
	Matrix<T> temporary(crsMatrix.rows(), crsMatrix.cols());
	for (SizeType i = 0; i < crsMatrix.rows(); i++) {
		//  for (SizeType k=0;k<crsMatrix.cols();k++) m(i,k)=0;
		for (int k = crsMatrix.getRowPtr(i); k < crsMatrix.getRowPtr(i + 1); k++)
			temporary(i, crsMatrix.getCol(k)) = crsMatrix.getValue(k);
	}

	m = std::move(temporary);
}

//! Transforms a full matrix into a Compressed-Row-Storage (CRS) Matrix
// Use the constructor if possible
template <typename T> void fullMatrixToCrsMatrix(CrsMatrix<T>& crsMatrix, const Matrix<T>& a)
{
	const T  zval      = 0.0;
	SizeType rows      = a.rows();
	SizeType cols      = a.cols();
	SizeType non_zeros = rows * cols;

	const bool use_push = true;

	if (use_push) {
		// ------------------------------
		// avoid filling array with zeros
		// ------------------------------
		crsMatrix.resize(rows, cols);
		crsMatrix.reserve(non_zeros);
	} else {
		crsMatrix.resize(rows, cols, non_zeros);
	};

	SizeType counter = 0;
	for (SizeType i = 0; i < rows; ++i) {
		crsMatrix.setRow(i, counter);
		for (SizeType j = 0; j < cols; ++j) {
			const T& val = a(i, j);
			if (val == zval)
				continue;

			if (use_push) {
				crsMatrix.pushValue(val);
				crsMatrix.pushCol(j);
			} else {
				crsMatrix.setValues(counter, val);
				crsMatrix.setCol(counter, j);
			};
			++counter;
		}
	}

	crsMatrix.setRow(rows, counter);
	crsMatrix.checkValidity();
}

/** If order==false then
                creates B such that
   B_{i1+j1*nout,i2+j2*nout)=A(j1,j2)\delta_{i1,i2} if order==true then creates
   B such that B_{i1+j1*na,i2+j2*na)=A(i1,i2)\delta_{j1,j2} where na=rank(A)
          */

template <typename T, typename VectorLikeType>
typename EnableIf<IsVectorLike<VectorLikeType>::True
                      && Loki::TypeTraits<typename VectorLikeType::value_type>::isFloat,
                  void>::Type
externalProduct(CrsMatrix<T>&                             B,
                const CrsMatrix<T>&                       A,
                SizeType                                  nout,
                const VectorLikeType&                     signs,
                bool                                      order,
                const PsimagLite::Vector<SizeType>::Type& permutationFull)
{
	if (A.rows() > 0)
		A.checkValidity();
	// -------------------------------------
	//  B = kron(eye, A)   if (is_A_fastest)
	//  B = kron(A, eye)   otherwise
	// -------------------------------------
	SizeType nrow_a   = A.rows();
	SizeType ncol_a   = A.cols();
	SizeType n        = nout;
	SizeType nrow_eye = n;
	SizeType ncol_eye = n;
	SizeType nnz_a    = A.nonZeros();

	SizeType nrow_b = n * nrow_a;
	SizeType ncol_b = n * ncol_a;
	SizeType nnz_b  = n * nnz_a;

	B.resize(nrow_b, ncol_b, nnz_b);

	bool is_a_fastest = order;

	if (nrow_a != ncol_a)
		throw RuntimeError("externalProduct: matrices must be square\n");

	// -----------------------
	// setup row pointers in B
	// Note: if (is_A_fastest)  then
	//          B( [ia,ie], [ja,je] ) = A(ia,ja) * eye(ie,je)
	//       else
	//          B( [ie,ia], [je,ja] ) = A(ia,ja) * eye(ie,je)
	//       endif
	//
	//  where [ia,ie] = ia + ie * nrow_A,   [ja,je] = ja + je * ncol_A
	//        [ie,ia] = ie + ia * nrow_eye, [je,ja] = je + ja * ncol_eye
	// -----------------------

	// -------------------------------------------------
	// calculate the number of nonzeros in each row of B
	// -------------------------------------------------
	std::vector<int> nnz_b_row(nrow_b);

	assert(nrow_A * nrow_eye <= permutationFull.size());

	for (SizeType ia = 0; ia < nrow_a; ia++) {

		SizeType nnz_row = A.getRowPtr(ia + 1) - A.getRowPtr(ia);

		for (SizeType ie = 0; ie < nrow_eye; ie++) {
			SizeType ib   = (is_a_fastest) ? permutationFull[ia + ie * nrow_a]
			                               : permutationFull[ie + ia * nrow_eye];
			nnz_b_row[ib] = nnz_row;
		};
	};

	// -------------------------------
	// setup row pointers in matrix B
	// -------------------------------
	std::vector<SizeType> b_rowptr(nrow_b);

	SizeType ip = 0;
	for (SizeType ib = 0; ib < nrow_b; ib++) {

		b_rowptr[ib] = ip;
		B.setRow(ib, ip);

		ip += nnz_b_row[ib];
	};
	assert(ip == nnz_B);
	B.setRow(nrow_b, nnz_b);

	// ---------------------------
	// copy entries into matrix B
	// ---------------------------

	// ----------------------------------------------
	// single pass over non-zero entries of matrix A
	// ----------------------------------------------
	for (SizeType ia = 0; ia < nrow_a; ia++) {
		for (int k = A.getRowPtr(ia); k < A.getRowPtr(ia + 1); k++) {

			// --------------------
			// entry aij = A(ia,ja)
			// --------------------
			SizeType ja  = A.getCol(k);
			T        aij = A.getValue(k);

			for (SizeType ie = 0; ie < nrow_eye; ie++) {
				SizeType je = ie;

				SizeType ib
				    = (is_a_fastest) ? ia + ie * nrow_a : ie + ia * nrow_eye;

				SizeType jb
				    = (is_a_fastest) ? ja + je * ncol_a : je + ja * ncol_eye;

				// --------------------
				// entry bij = B(ib,jb)
				// --------------------
				int alpha = ie;
				T   bij   = (is_a_fastest) ? aij : aij * signs[alpha];

				SizeType ip = b_rowptr[permutationFull[ib]];

				assert(jb < permutationFull.size());
				B.setCol(ip, permutationFull[jb]);
				B.setValues(ip, bij);

				++b_rowptr[permutationFull[ib]];
			};
		};
	};

	if (nrow_b != 0)
		B.checkValidity();
}

//-------

/** If order==false then
                creates C such that C_{i1+j1*nout,i2+j2*nout)=A(j1,j2)B_{i1,i2}
                if order==true then
                creates C such that C_{i1+j1*na,i2+j2*na)=A(i1,i2)B_{j1,j2}
                where na=rank(A) and nout = rank(B)
          */

template <typename T, typename VectorLikeType>
typename EnableIf<IsVectorLike<VectorLikeType>::True
                      && Loki::TypeTraits<typename VectorLikeType::value_type>::isFloat,
                  void>::Type
externalProduct(CrsMatrix<T>&                             C,
                const CrsMatrix<T>&                       A,
                const CrsMatrix<T>&                       B,
                const VectorLikeType&                     signs,
                bool                                      order,
                const PsimagLite::Vector<SizeType>::Type& permutationFull)
{
	const SizeType nfull = permutationFull.size();

	Vector<SizeType>::Type perm(nfull);
	for (SizeType i = 0; i < nfull; ++i)
		perm[permutationFull[i]] = i;

	const SizeType      nout       = B.rows();
	const SizeType      na         = A.rows();
	const SizeType      nout_or_na = (!order) ? nout : na;
	const CrsMatrix<T>& aor_b      = (!order) ? A : B;
	const CrsMatrix<T>& bor_a      = (!order) ? B : A;
	assert(A.rows() == A.cols());
	assert(B.rows() == B.cols());
	assert(nout * na == nfull);
	assert(signs.size() == noutOrNa);
	C.resize(nfull, nfull);
	SizeType counter = 0;
	for (SizeType i = 0; i < nfull; ++i) {
		C.setRow(i, counter);
		const SizeType ind = perm[i];
		ldiv_t         q   = std::ldiv(ind, nout_or_na);
		for (int k1 = bor_a.getRowPtr(q.rem); k1 < bor_a.getRowPtr(q.rem + 1); ++k1) {
			const SizeType col1 = bor_a.getCol(k1);
			for (int k2 = aor_b.getRowPtr(q.quot); k2 < aor_b.getRowPtr(q.quot + 1);
			     ++k2) {
				const SizeType col2 = aor_b.getCol(k2);
				SizeType       j    = permutationFull[col1 + col2 * nout_or_na];
				C.pushCol(j);
				C.pushValue(bor_a.getValue(k1) * aor_b.getValue(k2) * signs[q.rem]);
				++counter;
			}
		}
	}

	C.setRow(nfull, counter);
	C.checkValidity();
}

template <typename T>
void printFullMatrix(const CrsMatrix<T>& s,
                     const String&       name,
                     SizeType            how = 0,
                     double              eps = 1e-20)
{
	Matrix<T> fullm(s.rows(), s.cols());
	crsMatrixToFullMatrix(fullm, s);
	std::cout << "--------->   " << name;
	std::cout << " rank=" << s.rows() << "x" << s.cols() << " <----------\n";
	try {
		if (how == 1)
			mathematicaPrint(std::cout, fullm);
		if (how == 2)
			symbolicPrint(std::cout, fullm);
	} catch (std::exception& e) { }

	if (how == 0)
		fullm.print(std::cout, eps);
}

//! C = A*B,  all matrices are CRS matrices
template <typename S, typename S3, typename S2>
void multiply(CrsMatrix<S>& C, CrsMatrix<S3> const& A, CrsMatrix<S2> const& B)
{
	if (A.cols() != B.rows())
		throw RuntimeError("CrsMatrix::multiply: incompatible matrix dimensions\n");

	CrsMatrix<S> result;
	result.resize(A.rows(), B.cols());

	std::unordered_map<SizeType, S> accumulator;
	std::vector<SizeType>           columns;
	SizeType                        nonzeros = 0;

	for (SizeType row = 0; row < A.rows(); ++row) {
		result.setRow(row, nonzeros);
		accumulator.clear();
		columns.clear();

		for (int ka = A.getRowPtr(row); ka < A.getRowPtr(row + 1); ++ka) {
			const SizeType inner = A.getCol(ka);
			for (int kb = B.getRowPtr(inner); kb < B.getRowPtr(inner + 1); ++kb) {
				const SizeType column   = B.getCol(kb);
				const S        product  = A.getValue(ka) * B.getValue(kb);
				auto           inserted = accumulator.emplace(column, product);
				if (inserted.second)
					columns.push_back(column);
				else
					inserted.first->second += product;
			}
		}

		std::sort(columns.begin(), columns.end());
		for (const SizeType column : columns) {
			const S& value = accumulator.at(column);
			if (value == S(0))
				continue;

			result.pushCol(column);
			result.pushValue(value);
			++nonzeros;
		}
	}

	result.setRow(A.rows(), nonzeros);
	if (result.rows() > 0 && result.cols() > 0)
		result.checkValidity();
	C.swap(result);
}

// vector2 = sparseMatrix * vector1
template <class S>
void multiply(typename Vector<S>::Type&       v2,
              const CrsMatrix<S>&             m,
              const typename Vector<S>::Type& v1)
{
	SizeType n = m.rows();
	v2.resize(n);
	for (SizeType i = 0; i < n; i++) {
		v2[i] = 0;
		for (int j = m.getRowPtr(i); j < m.getRowPtr(i + 1); j++) {
			v2[i] += m.getValue(j) * v1[m.getCol(j)];
		}
	}
}

//! Sets B=transpose(conjugate(A))
template <typename S, typename S2> void transposeConjugate(CrsMatrix<S>& B, const CrsMatrix<S2>& A)
{
	SizeType nrow_a = A.rows();
	SizeType ncol_a = A.cols();
	SizeType nrow_b = ncol_a;
	SizeType ncol_b = nrow_a;

	SizeType nnz_a = A.nonZeros();
	SizeType nnz_b = nnz_a;

	B.resize(nrow_b, ncol_b, nnz_b);

	std::vector<SizeType> nnz_count(ncol_a, 0);

	// ----------------------------------------------------
	// 1st pass to count number of nonzeros per column in A
	// which is equivalent to number of nonzeros
	// per row in B = transpose(conjugate(A))
	// ----------------------------------------------------
	for (SizeType ia = 0; ia < nrow_a; ia++) {
		for (int k = A.getRowPtr(ia); k < A.getRowPtr(ia + 1); k++) {
			SizeType ja = A.getCol(k);
			++nnz_count[ja];
		};
	};

	// -----------------------
	// setup row pointers in B
	// -----------------------
	SizeType ipos = 0;
	for (SizeType ib = 0; ib < nrow_b; ib++) {
		B.setRow(ib, ipos);
		ipos += nnz_count[ib];
	};
	assert(ipos == nnz_B);
	B.setRow(nrow_b, nnz_b);

	// -------------------------------
	// setup row pointers in B matrix
	// -------------------------------
	std::vector<SizeType> b_rowptr(nrow_b);
	for (SizeType ib = 0; ib < nrow_b; ib++) {
		b_rowptr[ib] = B.getRowPtr(ib);
	};

	// -------------------------------------------
	// 2nd pass over matrix A to assign values to B
	// -------------------------------------------
	for (SizeType ia = 0; ia < nrow_a; ia++) {
		for (int k = A.getRowPtr(ia); k < A.getRowPtr(ia + 1); k++) {

			SizeType ja  = A.getCol(k);
			S2       aij = A.getValue(k);

			// ---------------------------------
			// B(ib=ja,jb=ia) = conj( A(ia,ja) )
			// ---------------------------------
			SizeType ib = ja;
			SizeType jb = ia;

			SizeType ip = b_rowptr[ib];

			// B.colind_[ ip ] = jb;
			// B.values_[ ip ] = PsimagLite::conj( aij );

			B.setCol(ip, jb);
			B.setValues(ip, PsimagLite::conj(aij));

			++b_rowptr[ib];
		};
	};
}

//! Sets A=B*b1+C*c1, restriction: B.size has to be larger or equal than C.size
template <typename T, typename T1>
void operatorPlus(CrsMatrix<T>& A, const CrsMatrix<T>& B, T1& b1, const CrsMatrix<T>& C, T1& c1)
{
	const T zero = static_cast<T>(0.0);

	SizeType nrow_b = B.rows();
	SizeType ncol_b = B.cols();
	SizeType nrow_c = C.rows();
	SizeType ncol_c = C.cols();

	// ------------------------------
	// nrow_A = std::max( nrow_B, nrow_C )
	// ncol_A = std::max( ncol_B, ncol_C )
	// ------------------------------
	SizeType nrow_a = (nrow_b >= nrow_c) ? nrow_b : nrow_c;
	SizeType ncol_a = (ncol_b >= ncol_c) ? ncol_b : ncol_c;

	A.resize(nrow_a, ncol_a);

	// ------------------------------------------------------
	// TODO: using A.resize(nrow_A,ncol_A,nnz_A) may not work correctly
	// ------------------------------------------------------
	const bool set_nonzeros = true;
	if (set_nonzeros) {
		SizeType nnz_b = B.nonZeros();
		SizeType nnz_c = C.nonZeros();

		// -----------------------------------------------
		// worst case when no overlap in sparsity pattern
		// between matrix B and matrix C
		// -----------------------------------------------
		SizeType nnz_a = nnz_b + nnz_c;

		A.reserve(nnz_a);
	};

	// ------------------------------------------
	// temporary vectors to accelerate processing
	// ------------------------------------------
	std::vector<T>    value_tmp(ncol_a, zero);
	std::vector<bool> is_examined_already(ncol_a, false);

	SizeType counter = 0;
	for (SizeType irow = 0; irow < nrow_a; irow++) {
		A.setRow(irow, counter);

		const bool is_valid_b_row = (irow < nrow_b);
		const bool is_valid_c_row = (irow < nrow_c);

		const int kstart_b = (is_valid_b_row) ? B.getRowPtr(irow) : 0;
		const int kend_b   = (is_valid_b_row) ? B.getRowPtr(irow + 1) : 0;

		const int kstart_c = (is_valid_c_row) ? C.getRowPtr(irow) : 0;
		const int kend_c   = (is_valid_c_row) ? C.getRowPtr(irow + 1) : 0;

		// --------------------------------
		// check whether there is work to do
		// --------------------------------
		const bool has_work = ((kend_b - kstart_b) + (kend_c - kstart_c) >= 1);
		if (!has_work)
			continue;

		// -------------------------------
		// add contributions from matrix B and matrix C
		// -------------------------------
		for (int k = kstart_b; k < kend_b; k++) {
			const T        bij  = B.getValue(k);
			const SizeType jcol = B.getCol(k);

			assert(jcol < ncol_A);

			value_tmp[jcol] += (bij * b1);
		};

		for (int k = kstart_c; k < kend_c; k++) {
			const T        cij  = C.getValue(k);
			const SizeType jcol = C.getCol(k);

			assert(jcol < ncol_A);

			value_tmp[jcol] += (cij * c1);
		};

		// --------------------
		// copy row to matrix A
		// --------------------

		for (int k = kstart_b; k < kend_b; k++) {
			const SizeType jcol = B.getCol(k);
			if (!is_examined_already[jcol]) {
				is_examined_already[jcol] = true;

				const T    aij     = value_tmp[jcol];
				const bool is_zero = (aij == zero);
				if (!is_zero) {
					A.pushCol(jcol);
					A.pushValue(aij);
					counter++;
				};
			};
		};

		for (int k = kstart_c; k < kend_c; k++) {
			const SizeType jcol = C.getCol(k);
			if (!is_examined_already[jcol]) {
				is_examined_already[jcol] = true;

				const T    aij     = value_tmp[jcol];
				const bool is_zero = (aij == zero);
				if (!is_zero) {
					A.pushCol(jcol);
					A.pushValue(aij);
					counter++;
				};
			};
		};

		// --------------------------------------------------
		// reset vectors valueTmp[] and is_examined_already[]
		// --------------------------------------------------

		for (int k = kstart_b; k < kend_b; k++) {
			const SizeType jcol       = B.getCol(k);
			value_tmp[jcol]           = zero;
			is_examined_already[jcol] = false;
		};

		for (int k = kstart_c; k < kend_c; k++) {
			const SizeType jcol       = C.getCol(k);
			value_tmp[jcol]           = zero;
			is_examined_already[jcol] = false;
		};

	}; // end for irow

	A.setRow(nrow_a, counter);

	// ----------------------------------------
	// set exact number of nonzeros in matrix A
	// ----------------------------------------
	SizeType nnz_a = counter;
	A.resize(nrow_a, ncol_a, nnz_a);

	A.checkValidity();
}

//! Sets A=B0*b0+B1*b1 + ...
template <typename T, typename T1>
void sum(CrsMatrix<T>&                           A,
         const std::vector<const CrsMatrix<T>*>& Bmats,
         const std::vector<T1>&                  bvec)
{
	SizeType bmats_size = Bmats.size();

	// ------------------------------
	// nrow_A = std::max( nrow_B(:) )
	// ncol_A = std::max( ncol_B(:) )
	// ------------------------------
	SizeType nrow_a  = 0;
	SizeType ncol_a  = 0;
	SizeType nnz_sum = 0;
	SizeType nnz_max = 0;

	for (SizeType imat = 0; imat < bmats_size; ++imat) {
		assert(imat < Bmats.size());
		const CrsMatrix<T>& this_mat = *(Bmats[imat]);
		SizeType            nrow_b   = this_mat.rows();
		SizeType            ncol_b   = this_mat.cols();
		SizeType            nnz_b    = this_mat.nonZeros();

		nrow_a = (nrow_b > nrow_a) ? nrow_b : nrow_a;
		ncol_a = (ncol_b > ncol_a) ? ncol_b : ncol_a;

		nnz_max = (nnz_b > nnz_max) ? nnz_b : nnz_max;
		nnz_sum += nnz_b;
	}

	A.resize(nrow_a, ncol_a);

	// ---------------------------------------------------
	// lower bound for total number of nonzeros is nnz_max
	// upper bound for total number of nonzeros is nnz_sum
	// initially set it to 2 * nnz_max
	// ---------------------------------------------------
	A.reserve((2 * nnz_max > nnz_sum) ? nnz_sum : 2 * nnz_max);

	// ------------------------------------------------------
	// TODO: using A.resize(nrow_A,ncol_A,nnz_A) may not work correctly
	// ------------------------------------------------------
	const bool set_nonzeros = true;
	if (set_nonzeros) {
		// -----------------------------------------------
		// worst case when no overlap in sparsity pattern
		// among matrices
		// -----------------------------------------------

		A.reserve(nnz_sum);
	}

	// ------------------------------------------
	// temporary vectors to accelerate processing
	// ------------------------------------------
	std::vector<T>    value_tmp(ncol_a);
	std::vector<bool> is_examined_already(ncol_a, false);

	std::vector<SizeType> column_index;
	column_index.reserve(ncol_a);

	SizeType counter = 0;
	for (SizeType irow = 0; irow < nrow_a; ++irow) {
		A.setRow(irow, counter);

		column_index.clear();

		for (SizeType imat = 0; imat < bmats_size; ++imat) {
			assert(imat < Bmats.size());
			const CrsMatrix<T>& this_mat       = *(Bmats[imat]);
			const SizeType      nrow_b         = this_mat.rows();
			const bool          is_valid_b_row = (irow < nrow_b);

			const SizeType kstart_b = (is_valid_b_row) ? this_mat.getRowPtr(irow) : 0;
			const SizeType kend_b = (is_valid_b_row) ? this_mat.getRowPtr(irow + 1) : 0;

			// --------------------------------
			// check whether there is work to do
			// --------------------------------
			const bool has_work = ((kend_b - kstart_b) >= 1);
			if (!has_work)
				continue;

			// -------------------------------
			// add contributions from matrix Bmats[i]
			// -------------------------------
			const T1 b1 = bvec[imat];
			for (SizeType k = kstart_b; k < kend_b; ++k) {
				const T        bij  = this_mat.getValue(k);
				const SizeType jcol = this_mat.getCol(k);

				assert(jcol < ncol_A);

				if (is_examined_already[jcol]) {
					value_tmp[jcol] += (bij * b1);
				} else {
					// ------------------------------------
					// new column entry not examined before
					// ------------------------------------

					is_examined_already[jcol] = true;

					value_tmp[jcol] = (bij * b1);

					column_index.push_back(jcol);
				}
			}
		} // end for imat

		// --------------------------------------
		// copy row to matrix A and reset vectors
		// --------------------------------------

		const SizeType kmax = column_index.size();
		for (SizeType k = 0; k < kmax; ++k) {
			const SizeType jcol = column_index[k];

			A.pushCol(jcol);
			A.pushValue(value_tmp[jcol]);
			is_examined_already[jcol] = false;
		}

		counter += kmax;
	} // end for irow

	SizeType nnz_a = counter;
	A.setRow(nrow_a, nnz_a);

	// ----------------------------------------
	// set exact number of nonzeros in matrix A
	//
	// note: this might be expensive in allocating another copy
	// and copying all the non-zero entries
	// ----------------------------------------
	const bool set_exact_nnz = true;
	if (set_exact_nnz)
		A.resize(nrow_a, ncol_a, nnz_a);

	A.checkValidity();
}

template <typename T> bool isHermitian(const CrsMatrix<T>& A, bool verbose = false)
{
	if (A.rows() != A.cols())
		return false;
	Matrix<T> dense;
	crsMatrixToFullMatrix(dense, A);
	return isHermitian(dense, verbose);
}

template <typename T> bool isAntiHermitian(const CrsMatrix<T>& A)
{
	if (A.rows() != A.cols())
		return false;
	Matrix<T> dense;
	crsMatrixToFullMatrix(dense, A);
	return isAntiHermitian(dense);
}

template <typename T>
void fromBlockToFull(CrsMatrix<T>& Bfull, const CrsMatrix<T>& B, SizeType offset)
{
	const bool use_push    = true;
	int        nrows_bfull = Bfull.rows();
	int        ncols_bfull = Bfull.cols();
	int        nnz_bfull   = B.nonZeros();
	Bfull.clear();

	if (use_push) {
		Bfull.resize(nrows_bfull, ncols_bfull);
		Bfull.reserve(nnz_bfull);
	} else {
		Bfull.resize(nrows_bfull, ncols_bfull, nnz_bfull);
	};

	int counter = 0;
	for (SizeType i = 0; i < offset; ++i)
		Bfull.setRow(i, counter);

	for (SizeType ii = 0; ii < B.rows(); ++ii) {
		SizeType i = ii + offset;
		Bfull.setRow(i, counter);
		for (int jj = B.getRowPtr(ii); jj < B.getRowPtr(ii + 1); ++jj) {
			SizeType j   = B.getCol(jj) + offset;
			T        tmp = B.getValue(jj);
			if (use_push) {
				Bfull.pushCol(j);
				Bfull.pushValue(tmp);
			} else {
				Bfull.setCol(counter, j);
				Bfull.setValues(counter, tmp);
			};
			counter++;
		}
	}

	for (SizeType i = B.rows() + offset; i < Bfull.rows(); ++i)
		Bfull.setRow(i, counter);

	Bfull.setRow(Bfull.rows(), counter);
	Bfull.checkValidity();
}

template <class T>
bool isDiagonal(const CrsMatrix<T>& A, double eps = 1e-6, bool checkForIdentity = false)
{
	if (A.rows() != A.cols())
		return false;
	SizeType n  = A.rows();
	const T  f1 = (-1.0);
	for (SizeType i = 0; i < n; i++) {
		for (int k = A.getRowPtr(i); k < A.getRowPtr(i + 1); k++) {
			SizeType col = A.getCol(k);
			const T& val = A.getValue(k);
			if (checkForIdentity && col == i && PsimagLite::norm(val + f1) > eps) {
				return false;
			}
			if (col != i && PsimagLite::norm(val) > eps) {
				return false;
			}
		}
	}
	return true;
}

template <class T> bool isTheIdentity(const CrsMatrix<T>& A, double eps = 1e-6)
{
	return isDiagonal(A, eps, true);
}

template <typename T> typename Real<T>::Type norm2(const CrsMatrix<T>& m)
{
	T val = 0;
	for (SizeType i = 0; i < m.values_.size(); i++)
		val += PsimagLite::conj(m.values_[i]) * m.values_[i];

	return PsimagLite::real(val);
}

template <typename T> Matrix<T> multiplyTc(const CrsMatrix<T>& a, const CrsMatrix<T>& b)
{

	CrsMatrix<T> bb, c;
	transposeConjugate(bb, b);
	multiply(c, a, bb);
	Matrix<T> cc;
	crsMatrixToFullMatrix(cc, c);
	return cc;
}

template <typename T> bool isZero(const CrsMatrix<T>& A, double eps = 0)
{
	SizeType n = A.values_.size();
	for (SizeType i = 0; i < n; ++i) {
		if (std::abs(A.values_[i]) > eps)
			return false;
	}

	return true;
}

} // namespace PsimagLite
/*@}*/
#endif
