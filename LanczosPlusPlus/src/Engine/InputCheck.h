// SPDX-FileCopyrightText: Copyright (c) 2009-2013, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [LanczosPlusPlus++, Version 1.0.0]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup LanczosPlusPlus */
/*@{*/

/*! \file InputCheck.h
 *
 *  InputChecking functions
 */
#ifndef LPP_INPUT_CHECK_H
#define LPP_INPUT_CHECK_H
#include <PsimagLite/Geometry/Geometry.h>
#include <PsimagLite/Options.h>
#include <stdexcept>
#include <vector>

namespace LanczosPlusPlus {

class InputCheck {

	typedef PsimagLite::Options::Readable OptionsReadableType;

public:

	InputCheck()
	    : optsReadable_(0)
	{ }

	~InputCheck()
	{
		if (optsReadable_ != 0)
			delete optsReadable_;
	}

	PsimagLite::String import() const
	{
		PsimagLite::String str = PsimagLite::Geometry<int, int, LanczosGlobals>::import();

		str += "vector hubbardU;\n";
		str += "vector potentialV;\n";
		str += "string! Model;\n";
		str += "string! SolverOptions;\n";
		str += "string! Version;\n";
		str += "integer! InfiniteLoopKeptStates;\n";
		str += "string! OutputFile;\n";
		str += "matrix.integer FiniteLoops;\n";
		str += "integer RepeatFiniteLoopsFrom;\n";
		str += "integer RepeatFiniteLoopsTimes;\n";
		str += "integer TargetElectronsUp;\n";
		str += "integer TargetElectronsDown;\n";
		str += "integer TargetElectronsTotal;\n";
		str += "real GsWeight;\n";
		str += "real TSPTau;\n";
		str += "integer TSPTimeSteps;\n";
		str += "integer TSPAdvanceEach;\n";
		str += "string TSPAlgorithm;\n";
		str += "vector.integer TSPSites;\n";
		str += "vector.integer TSPLoops;\n";
		str += "string TSPProductOrSum;\n";
		str += "string TSPOperator;\n";
		str += "string OperatorExpression;\n";
		str += "integer Threads = 1;\n";
		str += "integer Orbitals = 1;\n";
		str += "string FeAsMode;\n";
		str += "integer TargetSpinTimesTwo;\n";
		str += "integer UseSu2Symmetry;\n";
		str += "integer Pvectors;\n";
		str += "string TruncationTolerance;\n";
		str += "integer HeisenbergTwiceS;\n";
		str += "integer TargetSzPlusConst;\n";
		str += "integer SpinTwiceS;\n";
		str += "integer OrbitalTwiceS;\n";
		str += "real LambdaOne;\n";
		str += "real LambdaTwo;\n";
		str += "real CorrectionA;\n";
		str += "string RestartFilename;\n";
		str += "real LanczosEps;\n";
		str += "real TridiagonalEps;\n";
		str += "integer DynamicDmrgType;\n";
		str += "real CorrectionVectorFreqType;\n";
		str += "real CorrectionVectorEta;\n";
		str += "string CorrectionVectorAlgorithm;\n";
		str += "real CorrectionVectorOmega;\n";
		str += "string Intent;\n";
		str += "integer OpOnSiteThreshold;\n";
		str += "integer FirstRitz;\n";
		str += "integer CVnForFraction;\n";
		str += "real AnisotropyD;\n";
		str += "string FindSymmetrySector;\n";
		str += "string AddOnSiteHamiltonian;\n";
		str += "vector MagneticFieldX;\n";
		str += "vector MagneticFieldZ;\n";

		return str;
	}

	bool check(const PsimagLite::String&                           label,
	           const PsimagLite::Vector<PsimagLite::String>::Type& vec,
	           SizeType                                            line) const
	{
		if (label == "JMVALUES") {
			if (vec.size() != 2)
				return error1("JMVALUES", line);
			return true;
		} else if (label == "RAW_MATRIX" || label == "SpinOrbit") {
			SizeType row = atoi(vec[0].c_str());
			SizeType col = atoi(vec[1].c_str());
			SizeType n   = row * col;
			if (vec.size() != n + 2)
				return error1("RAW_MATRIX", line);
			return true;
		} else if (label == "Connectors") {
			return true;
		} else if (label == "MagneticField") {
			return true;
		} else if (label == "FiniteLoops") {
			SizeType n = atoi(vec[0].c_str());
			if (vec.size() != 3 * n + 1)
				return error1("FiniteLoops", line);
			return true;
		}
		return false;
	}

	bool checkSimpleLabel(const PsimagLite::String& label, SizeType line) const
	{
		// FIXME: needs implementation
		return true;
	}

	void check(const PsimagLite::String& label, const PsimagLite::String& val, SizeType)
	{
		if (label != "SolverOptions")
			return;
		PsimagLite::Vector<PsimagLite::String>::Type registerOpts;

		/* PSIDOC LanczosSolverOptions
		\begin{itemize}
		\item[none] Use this as a placeholder. ``none'' does not disable other options.
		\item[InternalProductStored] Stored the sparse matrix in memory before diagonalizing
		it.
		\item[InternalProductOnTheFly] Compute the sparse matrix on-the-fly while
		diagonalizing it.
		\item[printmatrix] Print the Hamiltonian matrix.
		\item[dumpmatrix] Use exact diagonalization instead of Lanczos diagonalization,
		and output all information to obtain the full spectrum.
		\item [setAffinities] TBW
		\end{itemize}
		*/
		registerOpts.push_back("none");
		registerOpts.push_back("InternalProductStored");
		registerOpts.push_back("InternalProductOnTheFly");
		registerOpts.push_back("printmatrix");
		registerOpts.push_back("dumpmatrix");
		registerOpts.push_back("setAffinities");

		PsimagLite::Options::Writeable optWriteable(
		    registerOpts, PsimagLite::Options::Writeable::PERMISSIVE);
		optsReadable_ = new OptionsReadableType(optWriteable, val);
	}

	bool isSet(const PsimagLite::String& thisOption) const
	{
		return optsReadable_->isSet(thisOption);
	}

	void usage(const char* progName)
	{
		std::cerr << "Usage: " << progName << " [-g -c] -f filename\n";
	}

private:

	bool error1(const PsimagLite::String& message, SizeType line) const
	{
		PsimagLite::String s(__FILE__);
		s += " : Input error for label " + message + " near line " + ttos(line) + "\n";
		throw PsimagLite::RuntimeError(s.c_str());
	}

	OptionsReadableType* optsReadable_;

}; // class InputCheck
} // namespace LanczosPlusPlus

/*@}*/
#endif
