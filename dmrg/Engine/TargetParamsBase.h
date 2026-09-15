// SPDX-FileCopyrightText: Copyright (c) 2014, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file TargetParamsBase.h
 *
 *  FIXME
 */
#ifndef TARGET_PARAMS_BASE_H
#define TARGET_PARAMS_BASE_H
#include <PsimagLite/FreqEnum.h>
#include <PsimagLite/MemResolv.h>
#include <stdexcept>
#include <vector>

namespace Dmrg {

template <typename ModelType> class TargetParamsBase {
public:

	using RealType           = typename ModelType::RealType;
	using PairFreqType       = std::pair<PsimagLite::FreqEnum, RealType>;
	using OperatorType       = typename ModelType::OperatorType;
	using PairType           = typename OperatorType::PairType;
	using SparseMatrixType   = typename OperatorType::SparseMatrixType;
	using ComplexOrRealType  = typename SparseMatrixType::value_type;
	using MatrixType         = PsimagLite::Matrix<ComplexOrRealType>;
	using VectorSizeType     = PsimagLite::Vector<SizeType>::Type;
	using VectorMatrixType   = typename PsimagLite::Vector<MatrixType>::Type;
	using VectorOperatorType = typename PsimagLite::Vector<OperatorType>::Type;
	using VectorRealType     = typename PsimagLite::Vector<RealType>::Type;

	enum class ConcatEnum
	{
		PRODUCT,
		SUM
	};

	enum class AlgorithmEnum
	{
		KRYLOV,
		CONJUGATE_GRADIENT,
		CHEBYSHEV,
		KRYLOVTIME,
		RUNGE_KUTTA,
		SUZUKI_TROTTER
	};

	TargetParamsBase(PsimagLite::String targeting)
	    : targeting_(targeting)
	{ }

	virtual ~TargetParamsBase() { }

	virtual SizeType sites() const = 0;

	virtual SizeType sectorIndex() const = 0;

	virtual SizeType levelIndex() const = 0;

	virtual SizeType sites(SizeType) const
	{
		PsimagLite::String s = "TargetParamsBase: unimplemented sites\n";
		throw PsimagLite::RuntimeError(s);
	}

	virtual void setOperator(SizeType, SizeType, const OperatorType&)
	{
		unimplemented("type(setOperator)");
	}

	virtual const VectorSizeType& startingLoops() const
	{
		PsimagLite::String s = "TargetParamsBase: unimplemented startingLoops\n";
		throw PsimagLite::RuntimeError(s);
	}

	virtual const VectorRealType& chebyTransform() const
	{
		PsimagLite::String s = "TargetParamsBase: unimplemented chebyTransform\n";
		throw PsimagLite::RuntimeError(s);
	}

	virtual ConcatEnum concatenation() const
	{
		throw PsimagLite::RuntimeError("concatenation");
	}

	virtual const VectorOperatorType& aOperators() const
	{
		PsimagLite::String s = "TargetParamsBase: unimplemented aOperators\n";
		throw PsimagLite::RuntimeError(s);
	}

	virtual RealType correctionA() const { return 0; }

	virtual SizeType type() const { return unimplementedInt("type"); }

	virtual void type(SizeType) { unimplemented("type(SizeType)"); }

	virtual SizeType advanceEach() const { return 0; }

	virtual SizeType cgSteps() const { return unimplementedInt("cgSteps"); }

	virtual PairFreqType omega() const
	{
		PsimagLite::String s("TargetParamsBase: unimplemented omega \n");
		throw PsimagLite::RuntimeError(s);
	}

	virtual void omega(PsimagLite::FreqEnum, RealType) { unimplemented("omega(RealType)"); }

	virtual RealType eta() const { return unimplemented("eta"); }

	virtual RealType cgEps() const { return unimplemented("cgEps"); }

	virtual AlgorithmEnum algorithm() const { throw PsimagLite::RuntimeError("algorithm"); }

	virtual RealType tau() const { return unimplemented("tau"); }

	virtual RealType maxTime() const { return unimplemented("maxTime"); }

	virtual VectorRealType& times()
	{
		throw PsimagLite::RuntimeError("times() unimplemented\n");
	}

	virtual const VectorRealType& times() const
	{
		throw PsimagLite::RuntimeError("times() unimplemented\n");
	}

	virtual bool noOperator() const { return static_cast<bool>(unimplemented("noOperator")); }

	virtual void noOperator(bool) { unimplemented("noOperator"); }

	virtual bool skipTimeZero() const
	{
		return static_cast<bool>(unimplemented("skipTimeZero"));
	}

	virtual bool isEnergyForExp() const
	{
		return static_cast<bool>(unimplemented("isEnergyForExp"));
	}

	virtual RealType energyForExp() const { return unimplemented("energyForExp"); }

	virtual RealType gsWeight() const { return unimplemented("gsWeight"); }

	virtual RealType timeDirection() const { return unimplemented("timeDirection"); }

	virtual PsimagLite::String targeting() const { return targeting_; }

private:

	RealType unimplemented(PsimagLite::String s) const
	{
		s = "TargetParamsBase: unimplemented " + s + "\n";
		throw PsimagLite::RuntimeError(s);
	}

	SizeType unimplementedInt(PsimagLite::String s) const
	{
		s = "TargetParamsBase: unimplemented " + s + "\n";
		throw PsimagLite::RuntimeError(s);
	}

	PsimagLite::String targeting_;
}; // class TargetParamsBase

} // namespace Dmrg

/*@}*/
#endif // TARGET_PARAMS_BASE_H
