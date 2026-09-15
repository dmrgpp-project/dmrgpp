// SPDX-FileCopyrightText: Copyright (c) 2009-2013, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file TridiagonalMatrix.h
 *
 *  A struct represent a tridiagonal matrix
 */
#ifndef TRIDIAGONAL_MATRIX_H
#define TRIDIAGONAL_MATRIX_H
#include "Matrix.h"

namespace PsimagLite {

template <typename FieldType> class TridiagonalMatrix {

	using RealType = typename Real<FieldType>::Type;

	static const bool diagWithLapack_ = true;

public:

	using VectorType     = typename Vector<FieldType>::Type;
	using VectorRealType = typename Vector<RealType>::Type;
	using value_type     = FieldType;

	TridiagonalMatrix() { }

	template <typename IoInputType> TridiagonalMatrix(IoInputType& io)
	{
		io.read(a_, "#Avector");
		io.read(b_, "#Bvector");
	}

	template <typename IoOutputType> void write(IoOutputType& io) const
	{
		io.write(a_, "#Avector");
		io.write(b_, "#Bvector");
	}

	void resize(SizeType n, FieldType value)
	{
		resize(n);
		for (SizeType i = 0; i < n; i++)
			a_[i] = b_[i] = value;
	}

	void resize(SizeType n)
	{
		assert(n > 0);
		a_.resize(n);
		b_.resize(n);
	}

	FieldType& a(SizeType i)
	{
		assert(i < a_.size());
		return a_[i];
	}

	FieldType& b(SizeType i)
	{
		assert(i < b_.size());
		return b_[i];
	}

	const FieldType& a(SizeType i) const
	{
		assert(i < a_.size());
		return a_[i];
	}

	const FieldType& b(SizeType i) const
	{
		assert(i < b_.size());
		return b_[i];
	}

	template <typename SomeMatrixType>
	void buildDenseMatrix(SomeMatrixType& m_dest, SizeType n = 0) const
	{
		if (n == 0)
			n = a_.size();
		SomeMatrixType m(n, n);
		for (SizeType i = 0; i < n - 1; ++i) {
			m(i, i)     = a_[i];
			m(i, i + 1) = b_[i + 1];
			m(i + 1, i) = PsimagLite::conj(b_[i + 1]);
		}

		m(n - 1, n - 1) = a_[n - 1];

		m_dest = std::move(m);
	}

	void diag(VectorRealType& eigs, SizeType nn) const
	{
		if (diagWithLapack_)
			diag2(eigs, nn);
		else
			ground(eigs, nn);
	}

	void push(const FieldType& a, const FieldType& b)
	{
		a_.push_back(a);
		b_.push_back(b);
	}

	SizeType size() const { return a_.size(); }

private:

	void diag2(VectorRealType&, SizeType) const;

	void ground(VectorRealType& groundD, SizeType nn) const
	{
		const int n = nn;

		groundD.resize(n);
		groundAllocations(n);

		const long int maxCounter = 10000;

		assert(a_.size() >= nn && b_.size() >= nn);
		for (SizeType i = 0; i < nn; ++i) {
			groundD[i]  = a_[i];
			groundE_[i] = b_[i];
		}

		RealType s          = 0;
		long int intCounter = 0;
		int      m          = 0;
		int      l          = 0;
		for (; l < n; l++) {
			do {
				intCounter++;
				if (intCounter > maxCounter) {
					std::cerr << "lanczos: ground: "
					             "premature exit ";
					std::cerr << "(may indicate an "
					             "internal error)\n";
					break;
				}

				for (m = l; m < n - 1; m++) {
					RealType dd = fabs(groundD[m]) + fabs(groundD[m + 1]);
					if ((fabs(groundE_[m]) + dd) == dd)
						break;
				}

				if (m != l) {
					RealType g
					    = (groundD[l + 1] - groundD[l]) / (2.0 * groundE_[l]);
					RealType r = sqrt(g * g + 1.0);
					g          = groundD[m] - groundD[l]
					    + groundE_[l] / (g + (g >= 0 ? fabs(r) : -fabs(r)));
					RealType p = 0.0;
					RealType c = 1.0;
					int      i = m - 1;
					for (s = 1.0; i >= l; i--) {
						RealType f      = s * groundE_[i];
						RealType h      = c * groundE_[i];
						groundE_[i + 1] = (r = sqrt(f * f + g * g));
						if (r == 0.0) {
							groundD[i + 1] -= p;
							groundE_[m] = 0.0;
							break;
						}

						s              = f / r;
						c              = g / r;
						g              = groundD[i + 1] - p;
						r              = (groundD[i] - g) * s + 2.0 * c * h;
						groundD[i + 1] = g + (p = s * r);
						g              = c * r - h;
					}

					if (r == 0.0 && i >= l)
						continue;
					groundD[l] -= p;
					groundE_[l] = g;
					groundE_[m] = 0.0;
				}
			} while (m != l);
		}

		std::sort(groundD.begin(), groundD.end());

		if (intCounter > maxCounter)
			throw RuntimeError(String(__FILE__) + "::ground(): internal error\n");
	}

	void groundAllocations(SizeType n) const
	{
		if (groundE_.size() != n) {
			groundE_.clear();
			groundE_.resize(n);
		}
	}

	VectorRealType         a_;
	VectorRealType         b_;
	mutable VectorRealType groundE_;
}; // class TridiagonalMatrix
} // namespace PsimagLite

/*@}*/
#endif
