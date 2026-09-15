// SPDX-FileCopyrightText: Copyright (c) 2009-2011, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.0.0]
// [by G.A., Oak Ridge National Laboratory]


/** \ingroup PsimagLite */
/*@{*/

/*! \file ParametersForSolver.h
 *
 * Parameters to be able to use LanczosSolver and ChebyshevSolver together
 *
 * Note: This is a struct: Add no functions except ctor/dtor!
 */

#ifndef PARAMETERS_FOR_SOLVER_H
#define PARAMETERS_FOR_SOLVER_H
#include "Matrix.h"
#include "ProgressIndicator.h"
#include "Random48.h"
#include "TridiagonalMatrix.h"
#include "TypeToString.h"
#include "Vector.h"

namespace PsimagLite {

template <typename RealType_> struct ParametersForSolver {

	using RealType = RealType_;

	static const SizeType MaxLanczosSteps = 1000000; // max number of internal Lanczos steps
	static const SizeType LanczosSteps    = 200; // max number of external Lanczos steps

	ParametersForSolver()
	    : steps(LanczosSteps)
	    , minSteps(4)
	    , tolerance(1e-12)
	    , stepsForEnergyConvergence(MaxLanczosSteps)
	    , eigsForStop(1)
	    , options("")
	    , oneOverA(0)
	    , b(0)
	    , Eg(0)
	    , weight(0)
	    , isign(0)
	    , lotaMemory(false)
	{ }

	template <typename IoInputType>
	ParametersForSolver(IoInputType& io, String prefix, int ind = -1)
	    : steps(LanczosSteps)
	    , minSteps(4)
	    , tolerance(1e-12)
	    , stepsForEnergyConvergence(MaxLanczosSteps)
	    , eigsForStop(1)
	    , options("none")
	    , oneOverA(0)
	    , b(0)
	    , Eg(0)
	    , weight(0)
	    , isign(0)
	    , lotaMemory(true)
	{

		rabbitHole(steps, prefix, ind, "Steps", io);

		rabbitHole(minSteps, prefix, ind, "MinSteps", io);

		rabbitHole(tolerance, prefix, ind, "Eps", io);

		rabbitHole(stepsForEnergyConvergence, prefix, ind, "StepsForEnergyConvergence", io);

		rabbitHole(eigsForStop, prefix, ind, "EigsForStop", io);

		rabbitHole(options, prefix, ind, "Options", io);

		rabbitHole(oneOverA, prefix, ind, "OneOverA", io);

		rabbitHole(b, prefix, ind, "B", io);

		rabbitHole(Eg, prefix, ind, "Energy", io);

		int x = 0;
		rabbitHole(x, prefix, ind, "NoSaveLanczosVectors", io);
		lotaMemory = (x > 0) ? 0 : 1;
	}

	template <typename T, typename IoInputType>
	static void rabbitHole(T& t, String prefix, int ind, String postfix, IoInputType& io)
	{
		if (ind >= 0) {
			bunnie(t, prefix, ind, postfix, io);
		} else {
			hare(t, prefix, postfix, io);
		}
	}

	template <typename T, typename IoInputType>
	static void hare(T& t, String prefix, String postfix, IoInputType& io)
	{
		// if prefix + postfix exists use it
		try {
			io.readline(t, prefix + postfix + "=");
			return;
		} catch (std::exception&) { }
	}

	template <typename T, typename IoInputType>
	static void bunnie(T& t, String prefix, SizeType ind, String postfix, IoInputType& io)
	{
		// if prefix + ind + postfix exists --> use it and return
		try {
			io.readline(t, prefix + ttos(ind) + postfix + "=");
			return;
		} catch (std::exception&) { }

		// if prefix + jnd + postfix exists with jnd < ind --> use the
		// largest jnd and return
		for (SizeType i = 0; i < ind; ++i) {
			const SizeType jnd = ind - i - 1;
			try {
				io.readline(t, prefix + ttos(jnd) + postfix + "=");
				return;
			} catch (std::exception&) { }
		}

		hare(t, prefix, postfix, io);
	}

	SizeType steps;
	SizeType minSteps;
	RealType tolerance;
	SizeType stepsForEnergyConvergence;
	SizeType eigsForStop;
	String   options;
	RealType oneOverA, b;
	RealType Eg;
	RealType weight;
	int      isign;
	bool     lotaMemory;
}; // class ParametersForSolver
} // namespace PsimagLite

/*@}*/
#endif // PARAMETERS_FOR_SOLVER_H
