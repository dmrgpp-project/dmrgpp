// SPDX-FileCopyrightText: Copyright (c) 2009-2013, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 2.0.0]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file TargetParamsTimeVectors.h
 *
 *  This is a structure to represent the parameters of the TimeStep Evolution
 *  algorithm. Don't add functions to this class because
 *  this class's data is all public
 */
#ifndef TARGET_PARAMS_TIME_VECTORS_H
#define TARGET_PARAMS_TIME_VECTORS_H
#include "TargetParamsCommon.h"

namespace Dmrg {
// Coordinates reading of TargetSTructure from input file
template <typename ModelType> class TargetParamsTimeVectors : public TargetParamsCommon<ModelType> {

public:

	using BaseType       = TargetParamsCommon<ModelType>;
	using RealType       = typename ModelType::RealType;
	using VectorRealType = typename BaseType::VectorRealType;

	template <typename IoInputter>
	TargetParamsTimeVectors(IoInputter&        io,
	                        PsimagLite::String targeting,
	                        const ModelType&   model)
	    : BaseType(io, targeting, model)
	    , advanceEach_(0)
	    , algorithm_(BaseType::AlgorithmEnum::KRYLOV)
	    , tau_(0)
	    , timeDirection_(1.0)
	{
		/*PSIDOC TargetParamsTimeVectors
		\item[TSPTau] [RealType], $\tau$ for the Krylov,
		see \cite{re:alvarez11} Section II.B and II.C.
		\item[TSPTimeSteps] [Integer]  $n_v$ as defined in
		\cite{re:alvarez11} Section II.B
		\item[TSPAdvanceEach] [Integer] Number of sites to sweep before
		advancing to the next time.
		\item[TSPAlgorithm] [String] Either
		\verb!Krylov! or \verb!RungeKutta! or \verb!SuzukiTrotter!\\
		Note that SuzukiTrotter is currently very experimental and unsupported.
		*/

		if (targeting == "TargetingExpression")
			return;

		io.readline(tau_, "TSPTau=");
		SizeType timeSteps = 0;
		io.readline(timeSteps, "TSPTimeSteps=");
		times_.resize(timeSteps);
		io.readline(advanceEach_, "TSPAdvanceEach=");
		PsimagLite::String s = "";

		io.readline(s, "TSPAlgorithm=");
		setAlgorithm(nullptr, s, &io);

		try {
			io.readline(timeDirection_, "TSPTimeFactor=");
		} catch (std::exception&) { }
	}

	VectorRealType& times() override { return times_; }

	const VectorRealType& times() const override { return times_; }

	SizeType advanceEach() const override { return advanceEach_; }

	typename BaseType::AlgorithmEnum algorithm() const override { return algorithm_; }

	RealType tau() const override { return tau_; }

	RealType timeDirection() const override { return timeDirection_; }

	const VectorRealType& chebyTransform() const override { return chebyTransform_; }

	template <typename IoInputter>
	void setAlgorithm(VectorRealType* chebyTransform, PsimagLite::String s, IoInputter* io)
	{
		if (io && chebyTransform)
			err("setAlgorithm: incorrect call (1)\n");

		if (!io && !chebyTransform)
			err("setAlgorithm: incorrect call (2)\n");

		if (s == "RungeKutta" || s == "rungeKutta" || s == "rungekutta") {
			algorithm_ = BaseType::AlgorithmEnum::RUNGE_KUTTA;
		} else if (s == "SuzukiTrotter" || s == "suzukiTrotter" || s == "suzukitrotter") {
			algorithm_ = BaseType::AlgorithmEnum::SUZUKI_TROTTER;
		} else if (s == "Chebyshev") {
			algorithm_ = BaseType::AlgorithmEnum::CHEBYSHEV;

			if (io)
				io->read(chebyTransform_, "ChebyshevTransform");
			else
				chebyTransform_ = *chebyTransform;
			if (chebyTransform_.size() != 2)
				err("ChebyshevTransform must be a vector of two real entries\n");
		} else if (s == "Krylov") {
			algorithm_ = BaseType::AlgorithmEnum::KRYLOV;
		} else {
			err("Unknown algorithm " + s + "\n");
		}
	}

private:

	VectorRealType                   times_;
	SizeType                         advanceEach_;
	typename BaseType::AlgorithmEnum algorithm_;
	RealType                         tau_;
	RealType                         timeDirection_;
	VectorRealType                   chebyTransform_;
}; // class TargetParamsTimeVectors

template <typename ModelType>
inline std::ostream& operator<<(std::ostream& os, const TargetParamsTimeVectors<ModelType>& t)
{
	os << "TargetParams.type=TimeVectors";
	os << "TargetParams.tau=" << t.tau() << "\n";
	os << "TargetParams.timeSteps=" << t.timeSteps() << "\n";
	os << "TargetParams.advanceEach=" << t.advanceEach() << "\n";
	os << "TargetParams.algorithm=" << t.algorithm() << "\n";
	os << "TargetParams.timeDirection=" << t.timeDirection() << "\n";
	return os;
}
} // namespace Dmrg

/*@}*/
#endif // TARGET_PARAMS_TIME_VECTORS_H
