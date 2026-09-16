// SPDX-FileCopyrightText: Copyright (c) 2009-2013, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup PsimagLite */
/*@{*/

/*! \file ChebyshevSolver.h
 *
 *  A class to represent a generic Lanczos Solver
 *
 */

#ifndef CHEBYSHEV_SOLVER_H_
#define CHEBYSHEV_SOLVER_H_
#include "ChebyshevSerializer.h"
#include "LanczosSolver.h"
#include "Matrix.h"
#include "MatrixSolverBase.hh"
#include "ParametersForSolver.h"
#include "ProgressIndicator.h"
#include "Random48.h"
#include "TridiagonalMatrix.h"
#include "TypeToString.h"
#include "Vector.h"

namespace PsimagLite {

/** MatrixType must have the following interface:
 * RealType type to indicate the matrix type
 * rows() member function to indicate the rank of the matrix
 * matrixVectorProduct(typename Vector< RealType>::Type& x,const
 *        typename Vector< RealType>::Type& const y)
 *    member function that implements the operation x += Hy
 *
 */
template <typename MatrixType_> class ChebyshevSolver {

	using ComplexOrRealType   = typename MatrixType_::value_type;
	using RealType            = typename Real<ComplexOrRealType>::Type;
	using VectorType          = std::vector<ComplexOrRealType>;
	using LanczosVectorsType  = LanczosVectors<MatrixType_, VectorType>;
	using DenseMatrixType     = typename LanczosVectorsType::DenseMatrixType;
	using DenseMatrixRealType = typename LanczosVectorsType::DenseMatrixRealType;
	using VectorVectorType    = typename LanczosVectorsType::VectorVectorType;

public:

	using ParametersSolverType  = ParametersForSolver<RealType>;
	using MatrixType            = MatrixType_;
	using TridiagonalMatrixType = TridiagonalMatrix<RealType>;
	using PostProcType          = ChebyshevSerializer<TridiagonalMatrixType>;
	using RngType               = PsimagLite::Random48<RealType>;

	enum
	{
		WITH_INFO   = 1,
		DEBUG       = 2,
		ALLOWS_ZERO = 4
	};

	ChebyshevSolver(MatrixType const& mat, ParametersSolverType& params)
	    : progress_("ChebyshevSolver")
	    , mat_(mat)
	    , params_(params)
	    , mode_(WITH_INFO)
	    , rng_(343311)
	    , lanczosVectors_(
	          mat,
	          params.lotaMemory,
	          params.steps,
	          MatrixSolverBase<MatrixType>::isReorthoEnabled(params.options, params.lotaMemory))
	{
		params.steps = 400;
		setMode(params.options);
		computeAandB();
		PsimagLite::OstringStream                     msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg << "Constructing... mat.rank=" << mat_.rows() << " steps=" << params.steps;
		progress_.printline(msgg, std::cout);
	}

	void computeGroundState(RealType&, VectorType&) { unimplemented("computeGroundState"); }

	void computeGroundState(RealType& gsEnergy, VectorType& z, const VectorType& initialVector)
	{
		if (mode_ & DEBUG) {
			computeGroundStateTest(gsEnergy, z, initialVector);
			return;
		}
		unimplemented("computeGroundState");
	}

	void buildDenseMatrix(DenseMatrixType&, const TridiagonalMatrixType&) const
	{
		unimplemented("buildDenseMatrix");
	}

	void push(TridiagonalMatrixType& ab, const RealType& a, const RealType& b) const
	{
		ab.push_back(a);
		ab.push_back(b);
	}

	//! ab.a contains the even moments
	//! ab.b contains the odd moments
	void decomposition(const VectorType& initVector, TridiagonalMatrixType& ab)
	{
		VectorType x(initVector.size(), 0.0);
		VectorType y = initVector;

		lanczosVectors_.prepareMemory(y.size(), lanczosVectors_.cols());
		ab.resize(2 * params_.steps, 0);
		SizeType cols = lanczosVectors_.cols();
		for (SizeType j = 0; j < cols; ++j) {
			if (lanczosVectors_.lotaMemory())
				lanczosVectors_.saveVector(y, j);

			RealType atmp = 0;
			RealType btmp = 0;
			oneStepDec(x, y, atmp, btmp, j);
			ab.a(j) = 2 * atmp - ab.a(0);
			ab.b(j) = 2 * btmp - ab.b(0);
		}

		// lanczosVectors_.resize(cols); <--- not needed because all
		// steps are performed
		//                                    and there is no early exit
		//                                    here
	}

	//! atmp = < phi_n | phi_n>
	//! btmp = < phi_n | phi_{n+1}>
	void oneStepDec(VectorType& x,
	                VectorType& y,
	                RealType&   atmp,
	                RealType&   btmp,
	                SizeType    jind) const
	{
		bool       isFirst = (jind == 0);
		VectorType z(x.size(), 0.0);
		mat_.matrixVectorProduct(z, y); // z+= Hy
		// scale matrix:
		z -= params_.b * y;
		z *= params_.oneOverA;

		RealType val = (isFirst) ? 1.0 : 2.0;

		atmp = 0.0;
		for (SizeType i = 0; i < mat_.rows(); i++)
			atmp += PsimagLite::real(y[i] * PsimagLite::conj(y[i]));

		for (SizeType i = 0; i < mat_.rows(); i++) {
			ComplexOrRealType tmp = val * z[i] - x[i];
			x[i]                  = y[i];
			y[i]                  = tmp;
		}

		btmp = 0.0;
		for (SizeType i = 0; i < mat_.rows(); i++)
			btmp += PsimagLite::real(y[i] * PsimagLite::conj(x[i]));
	}

	SizeType steps() const { return params_.steps; }

	void lanczosVectorsSwap(DenseMatrixType& V)
	{
		DenseMatrixType* ptr = lanczosVectors_.data();
		if (!ptr)
			err("LanczosSolver::lanczosVectors() called but no "
			    "data stored\n");
		return ptr->swap(V);
	}

private:

	void unimplemented(const String& s) const
	{
		String s2("Hmmm...this ain't looking good...");
		s2 += String(__FILE__) + " " + ttos(__LINE__) + " ";
		s2 += s;
		throw RuntimeError(s);
	}

	void setMode(const String& options)
	{
		if (options.find("lanczosdebug") != String::npos)
			mode_ |= DEBUG;

		if (options.find("lanczosAllowsZero") != String::npos)
			mode_ |= ALLOWS_ZERO;
	}

	void info(RealType energyTmp, const VectorType& x, std::ostream& os)
	{
		RealType norma = norm(x);

		if (norma < 1e-5 || norma > 100) {
			std::cerr << "norma=" << norma << "\n";
			// throw RuntimeError("Norm\n");
		}

		PsimagLite::OstringStream                     msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg << "Found Energy=" << energyTmp << " after " << params_.steps;
		msg << " iterations, "
		    << " orig. norm=" << norma;
		progress_.printline(msgg, os);
	}

	//! only for debugging:
	void
	computeGroundStateTest(RealType& gsEnergy, VectorType& z, const VectorType& initialVector)
	{
		unimplemented("computeGroundStateTest");
	}

	class InternalMatrix {
	public:

		using value_type = ComplexOrRealType;

		InternalMatrix(const MatrixType& mat)
		    : matx_(mat)
		    , y_(matx_.rows())
		{ }

		SizeType rows() const { return matx_.rows(); }

		void matrixVectorProduct(VectorType& x, const VectorType& y) const
		{
			for (SizeType i = 0; i < y_.size(); i++)
				y_[i] = -y[i];
			matx_.matrixVectorProduct(x, y_);
		}

		ComplexOrRealType operator()(SizeType i, SizeType j) const { return matx_(i, j); }

	private:

		const MatrixType&  matx_;
		mutable VectorType y_;
	}; // class InternalMatrix

	void computeAandB()
	{
		PsimagLite::OstringStream                     msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg << "Asking LanczosSolver to compute spectrum bounds...";
		progress_.printline(msgg, std::cout);

		ParametersSolverType          params;
		InternalMatrix                mat2(mat_);
		RealType                      eMax = 0;
		LanczosSolver<InternalMatrix> lanczosSolver2(mat2, params);

		VectorType z2(mat_.rows(), 0);
		VectorType init(z2.size());
		PsimagLite::fillRandom(init);
		lanczosSolver2.computeOneState(eMax, z2, init, 0);

		VectorType                z(mat_.rows(), 0);
		LanczosSolver<MatrixType> lanczosSolver(mat_, params);
		RealType                  eMin = 0;
		lanczosSolver.computeOneState(eMin, z, init, 0);

		eMax = -eMax;
		eMax *= 3;
		eMin *= 3;
		assert(eMax - eMin > 1e-2);

		params_.oneOverA = 2.0 / (eMax - eMin);
		params_.b        = (eMax + eMin) / 2;

		PsimagLite::OstringStream                     msgg2(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg2 = msgg2();
		msg2 << "Spectrum bounds computed, eMax=" << eMax << " eMin=" << eMin;
		progress_.printline(msgg2, std::cout);
	}

	ProgressIndicator     progress_;
	MatrixType const&     mat_;
	ParametersSolverType& params_;
	SizeType              mode_;
	RngType               rng_;
	LanczosVectorsType    lanczosVectors_;
	//! Scaling factors for the Chebyshev expansion
}; // class ChebyshevSolver
} // namespace PsimagLite
/*@}*/
#endif // CHEBYSHEV_SOLVER_H_
