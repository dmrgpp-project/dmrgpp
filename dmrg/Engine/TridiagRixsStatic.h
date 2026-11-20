// SPDX-FileCopyrightText: Copyright (c) 2009,2017 UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 4.]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/
/** \file TridiagRixsStatic.h
 * This class performs the tridiagonalization of a modified Hamiltonian, which is
 * H0 + A_i for i the current site and A a local operator
 * This is needed to implement RixsStatic
 */
#ifndef TRIDIAGRIXSSTATIC_H
#define TRIDIAGRIXSSTATIC_H
#include "ApplyOperatorLocal.h"
#include <PsimagLite/LanczosSolver.h>

namespace Dmrg {

template <typename ModelType, typename LanczosSolverType, typename VectorWithOffsetType>
class TridiagRixsStatic {

	using ModelHelperType        = typename ModelType::ModelHelperType;
	using LeftRightSuperType     = typename ModelHelperType::LeftRightSuperType;
	using BasisWithOperatorsType = typename LeftRightSuperType::BasisWithOperatorsType;
	using OperatorsType          = typename BasisWithOperatorsType::OperatorsType;
	using OperatorType           = typename OperatorsType::OperatorType;
	using SparseMatrixType       = typename BasisWithOperatorsType::SparseMatrixType;
	using ComplexOrRealType      = typename SparseMatrixType::value_type;
	using RealType               = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using InputValidatorType     = typename ModelType::InputValidatorType;
	using ConcurrencyType        = PsimagLite::Concurrency;
	using SparseElementType      = typename SparseMatrixType::value_type;
	using MatrixType             = PsimagLite::Matrix<SparseElementType>;
	using ParametersSolverType   = typename LanczosSolverType::ParametersSolverType;
	using VectorWithOffsetType2  = VectorWithOffsets<ComplexOrRealType>;
	using ApplyOperatorLocalType
	    = ApplyOperatorLocal<LeftRightSuperType, VectorWithOffsetType2>;
	using VectorSizeType = typename VectorWithOffsetType2::VectorSizeType;

	class MyMatrixVector : public LanczosSolverType::LanczosMatrixType {

		using BasisType = typename LanczosSolverType::LanczosMatrixType;

	public:

		MyMatrixVector(ModelType const*                            model,
		               ModelHelperType const*                      modelHelper,
		               const OperatorType&                         A,
		               ProgramGlobals::DirectionEnum               dir,
		               typename ApplyOperatorLocalType::BorderEnum corner,
		               const VectorSizeType&                       weights)
		    : BasisType(model, modelHelper)
		    , applyOperatorLocal_(modelHelper->leftRightSuper())
		    , A_(A)
		    , dir_(dir)
		    , corner_(corner)
		    , fs_(modelHelper->leftRightSuper().left().electronsVector())
		    , // FIXME CHECK
		    x2_(weights, modelHelper->leftRightSuper().super())
		    , y2_(weights, modelHelper->leftRightSuper().super())
		{ }

		template <typename SomeVectorType>
		void matrixVectorProduct(SomeVectorType& x, SomeVectorType const& y) const
		{
			BasisType::matrixVectorProduct(x, y);
			MatrixType fullA;
			crsMatrixToFullMatrix(fullA, A_.data);
			if (isZero(fullA))
				return;
			// add here x += Ay
			x2_.setDataInSector(x, 0);
			y2_.setDataInSector(y, 0);
			applyOperatorLocal_(x2_, y2_, A_, fs_, dir_, corner_);
			x2_.extract(x, 0);
		}

	private:

		ApplyOperatorLocalType                      applyOperatorLocal_;
		const OperatorType&                         A_;
		ProgramGlobals::DirectionEnum               dir_;
		typename ApplyOperatorLocalType::BorderEnum corner_;
		FermionSign                                 fs_;
		mutable VectorWithOffsetType2               x2_;
		mutable VectorWithOffsetType2               y2_;
	}; // class MyMatrixVector

	using MyMatrixVectorType      = MyMatrixVector;
	using MyLanczosSolverType     = PsimagLite::LanczosSolver<MyMatrixVectorType>;
	using MyTridiagonalMatrixType = typename MyLanczosSolverType::TridiagonalMatrixType;

public:

	using TargetVectorType        = typename PsimagLite::Vector<ComplexOrRealType>::Type;
	using MatrixComplexOrRealType = PsimagLite::Matrix<ComplexOrRealType>;
	using VectorMatrixFieldType   = typename PsimagLite::Vector<MatrixComplexOrRealType>::Type;

	TridiagRixsStatic(const LeftRightSuperType&     lrs,
	                  const ModelType&              model,
	                  InputValidatorType&           io,
	                  SizeType                      site,
	                  ProgramGlobals::DirectionEnum direction)
	    : lrs_(lrs)
	    , model_(model)
	    , io_(io)
	    , A_(io, model, false, "RS:")
	    , direction_(direction)
	{
		SizeType numberOfSites = model.geometry().numberOfSites();

		int site2 = ProgramGlobals::findBorderSiteFrom(site, direction, numberOfSites);
		corner_   = (site2 >= 0) ? ApplyOperatorLocalType::BORDER_YES
		                         : ApplyOperatorLocalType::BORDER_NO;
	}

	void operator()(const VectorWithOffsetType& phi,
	                VectorMatrixFieldType&      T,
	                VectorMatrixFieldType&      V,
	                VectorSizeType&             steps)
	{
		for (SizeType ii = 0; ii < phi.sectors(); ++ii) {
			SizeType i = phi.sector(ii);
			steps[ii]  = triDiag(phi, T[ii], V[ii], i);
		}
	}

private:

	SizeType triDiag(const VectorWithOffsetType& phi,
	                 MatrixComplexOrRealType&    T,
	                 MatrixComplexOrRealType&    V,
	                 SizeType                    i0)
	{
		VectorSizeType weights(lrs_.super().partition(), 0);
		weights[i0]          = phi.effectiveSize(i0);
		SizeType p           = lrs_.super().findPartitionNumber(phi.offset(i0));
		SizeType threadNum   = 0;
		SizeType currentTime = 0;
		typename ModelType::ModelHelperType modelHelper(p, lrs_, currentTime, threadNum);
		typename MyLanczosSolverType::LanczosMatrixType lanczosHelper(
		    &model_, &modelHelper, A_, direction_, corner_, weights);

		ParametersSolverType params(io_, "Tridiag");
		params.lotaMemory = true;
		params.threadId   = threadNum;

		MyLanczosSolverType lanczosSolver(lanczosHelper, params, &V);

		MyTridiagonalMatrixType ab;
		SizeType                total = phi.effectiveSize(i0);
		TargetVectorType        phi2(total);
		phi.extract(phi2, i0);
		lanczosSolver.decomposition(phi2, ab);
		lanczosSolver.buildDenseMatrix(T, ab);
		return lanczosSolver.steps();
	}

	const LeftRightSuperType&                   lrs_;
	const ModelType&                            model_;
	InputValidatorType&                         io_;
	OperatorType                                A_;
	typename ApplyOperatorLocalType::BorderEnum corner_;
	ProgramGlobals::DirectionEnum               direction_;
}; // class TridiagRixsStatic
} // namespace Dmrg
#endif // TRIDIAGRIXSSTATIC_H
