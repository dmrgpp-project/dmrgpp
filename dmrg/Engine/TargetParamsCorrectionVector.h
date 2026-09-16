// SPDX-FileCopyrightText: Copyright (c) 2009-2016-2018, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 4.]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file TargetParamsCorrectionVector.h
 *
 *  This is a structure to represent the parameters of the
 *  Correction vector DMRG algorithm.
 *  Don't add functions to this class because
 *  this class's data is all public
 */
#ifndef TARGET_PARAMS_CORRECTION_V_H
#define TARGET_PARAMS_CORRECTION_V_H

#include "TargetParamsCommon.h"
#include <PsimagLite/FreqEnum.h>

namespace Dmrg {
// Coordinates reading of TargetSTructure from input file
template <typename ModelType>
class TargetParamsCorrectionVector : public TargetParamsCommon<ModelType> {

public:

	using BaseType         = TargetParamsCommon<ModelType>;
	using RealType         = typename ModelType::RealType;
	using PairFreqType     = typename BaseType::BaseType::PairFreqType;
	using OperatorType     = typename ModelType::OperatorType;
	using PairType         = typename OperatorType::PairType;
	using SparseMatrixType = typename OperatorType::StorageType;
	using ComplexOrReal    = typename SparseMatrixType::value_type;
	using MatrixType       = PsimagLite::Matrix<ComplexOrReal>;

	template <typename IoInputter>
	TargetParamsCorrectionVector(IoInputter&        io,
	                             PsimagLite::String targeting,
	                             const ModelType&   model)
	    : BaseType(io, targeting, model)
	    , cgSteps_(1000)
	    , firstRitz_(0)
	    , nForFraction_(1)
	    , advanceEach_(0)
	    , cgEps_(1e-6)
	{
		io.readline(correctionA_, "CorrectionA=");
		io.readline(type_, "DynamicDmrgType=");
		PsimagLite::String tmp;
		io.readline(tmp, "CorrectionVectorFreqType=");
		PsimagLite::FreqEnum freqEnum = PsimagLite::FreqEnum::REAL;

		if (tmp == "Matsubara") {
			freqEnum = PsimagLite::FreqEnum::MATSUBARA;
		} else if (tmp != "Real") {
			PsimagLite::String msg("CorrectionVectorFreqType");
			throw PsimagLite::RuntimeError(msg += "must be either Real or Matsubara\n");
		}

		RealType omega;
		io.readline(omega, "CorrectionVectorOmega=");
		omega_ = PairFreqType(freqEnum, omega);
		io.readline(eta_, "CorrectionVectorEta=");

		io.readline(tmp, "CorrectionVectorAlgorithm=");
		if (tmp == "Krylov") {
			algorithm_ = BaseType::AlgorithmEnum::KRYLOV;
		} else if (tmp == "ConjugateGradient") {
			algorithm_ = BaseType::AlgorithmEnum::CONJUGATE_GRADIENT;
		} else if (tmp == "Chebyshev") {
			algorithm_ = BaseType::AlgorithmEnum::CHEBYSHEV;
		} else if (tmp == "KrylovTime") {
			algorithm_ = BaseType::AlgorithmEnum::KRYLOVTIME;
		} else {
			PsimagLite::String str("TargetParamsCorrectionVector ");
			str += "Unknown algorithm " + tmp + "\n";
			throw PsimagLite::RuntimeError(str);
		}

		try {
			io.readline(cgSteps_, "ConjugateGradientSteps=");
		} catch (std::exception&) { }

		try {
			io.readline(cgEps_, "ConjugateGradientEps=");
		} catch (std::exception&) { }

		try {
			int x = 0;
			io.readline(x, "TSPUseQns=");
			err("TSPUseQns= is no longer needed, please delete it from the input "
			    "file\n");
		} catch (std::exception&) { }

		try {
			io.readline(firstRitz_, "FirstRitz=");
		} catch (std::exception&) { }

		try {
			io.readline(nForFraction_, "CVnForFraction=");
		} catch (std::exception&) { }

		if (nForFraction_ > 1)
			io.readline(advanceEach_, "TSPAdvanceEach=");

		if (freqEnum == PsimagLite::FreqEnum::MATSUBARA && firstRitz_ != 0)
			err("FirstRitz must be 0 for Matsubara\n");
	}

	RealType correctionA() const override { return correctionA_; }

	SizeType type() const override { return type_; }

	void type(SizeType x) override { type_ = x; }

	SizeType cgSteps() const override { return cgSteps_; }

	PairFreqType omega() const override { return omega_; }

	void omega(PsimagLite::FreqEnum freqEnum, RealType x) override
	{
		omega_ = PairFreqType(freqEnum, x);
	}

	RealType eta() const override { return eta_; }

	RealType cgEps() const override { return cgEps_; }

	typename BaseType::AlgorithmEnum algorithm() const override { return algorithm_; }

	virtual SizeType firstRitz() const { return firstRitz_; }

	virtual SizeType nForFraction() const { return nForFraction_; }

	SizeType advanceEach() const override { return advanceEach_; }

private:

	SizeType                         type_;
	typename BaseType::AlgorithmEnum algorithm_;
	SizeType                         cgSteps_;
	SizeType                         firstRitz_;
	SizeType                         nForFraction_;
	SizeType                         advanceEach_;
	RealType                         correctionA_;
	PairFreqType                     omega_;
	RealType                         eta_;
	RealType                         cgEps_;
}; // class TargetParamsCorrectionVector

template <typename ModelType>
inline std::ostream& operator<<(std::ostream& os, const TargetParamsCorrectionVector<ModelType>& t)
{
	os << "TargetParams.type=AdaptiveDynamic\n";
	const TargetParamsCommon<ModelType>& tp = t;
	os << tp;
	os << "DynamicDmrgType=" << t.type() << "\n";
	os << "CorrectionVectorOmega=" << t.omega() << "\n";
	os << "CorrectionVectorEta=" << t.eta() << "\n";
	os << "ConjugateGradientSteps" << t.cgSteps() << "\n";
	os << "ConjugateGradientEps" << t.cgEps() << "\n";
	os << "firstRitz" << t.firstRitz() << "\n";
	os << "nForFraction" << t.nForFraction() << "\n";

	return os;
}
} // namespace Dmrg

/*@}*/
#endif // TARGET_PARAMS_CORRECTION_V_H
