// SPDX-FileCopyrightText: Copyright (c) 2009-2015, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup DMRG */
/*@{*/

/*! \file ToolBox.h
 *
 *
 */

#ifndef DMRG_TOOLBOX_H
#define DMRG_TOOLBOX_H
#include "ProgramGlobals.h"
#include <PsimagLite/Geometry/Geometry.h>
#include <PsimagLite/PsimagLite.h>
#include <PsimagLite/Vector.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

namespace Dmrg {

template <typename DmrgParametersType, typename GeometryType> class ToolBox {

	using PairSizeStringType = std::pair<SizeType, PsimagLite::String>;

	class GrepForLabel {

		using LongType = long int;

		struct InternalName {
			InternalName(PsimagLite::String label_, bool cooked_)
			    : cooked(cooked_)
			    , label(label_)
			{ }

			bool               cooked;
			PsimagLite::String label;
		}; // struct InternalName

	public:

		using ParametersType = InternalName;

		static void hook(std::ifstream& fin,
		                 PsimagLite::String,
		                 LongType              len,
		                 const ParametersType& params)
		{
			LongType          len2      = len;
			LongType          bufferLen = 1;
			std::stringstream ss;
			char*             buffer = new char[bufferLen];
			while (len2 >= bufferLen && !fin.eof()) {
				fin.read(buffer, bufferLen);
				ss << buffer[0];
				if (buffer[0] == '\n') {
					procLine(ss.str(), params);
					ss.str("");
				}

				if (len > 1)
					len2 -= bufferLen;
			}

			delete[] buffer;
		}

	private:

		static void procLine(PsimagLite::String line, const ParametersType& params)
		{
			if (line.find(params.label) == PsimagLite::String::npos)
				return;
			if (params.cooked)
				cookThisLine(line);
			else
				std::cout << line;
		}

		static void cookThisLine(PsimagLite::String line)
		{
			size_t             index = line.find(" after");
			PsimagLite::String line2 = line;
			if (index != PsimagLite::String::npos)
				line2 = line.erase(index, line.length());
			line                     = line2;
			PsimagLite::String magic = "eigenvalue= ";
			index                    = line.find(magic);
			if (index != PsimagLite::String::npos)
				line.erase(0, index + magic.length());
			std::cout << line << "\n";
		}
	}; // GrepForLabel

public:

	enum ActionEnum
	{
		ACTION_UNKNOWN,
		ACTION_GREP,
		ACTION_FILES,
		ACTION_INPUT,
		ACTION_ANALYSIS
	};

	using ParametersForGrepType = typename GrepForLabel::ParametersType;

	/* PSIDOC ToolBoxActions
	 Actions marked with an asterix are only meaningful postprocessing.
	 \begin{itemize}
	 \item[energy] or Energy or energies or Energies. It lists energies of all stages. (*)
	 \item[files] TBW
	 \item[input] It echoes the input file.
	 \item[analysis] or analyze. It opines about the needed ``m'' values for this run,
	 as well as the needed CPU and RAM that will be required.
	 \end{itemize}
	 */
	static ActionEnum actionCanonical(PsimagLite::String action)
	{
		if (action == "energy" || action == "Energy" || action == "energies"
		    || action == "Energies" || action == "grep")
			return ACTION_GREP;
		if (action == "files")
			return ACTION_FILES;
		if (action == "input")
			return ACTION_INPUT;
		if (action == "analysis" || action == "analyze")
			return ACTION_ANALYSIS;
		return ACTION_UNKNOWN;
	}

	static PsimagLite::String actions() { return "energies | grep | files | input |analysis"; }

	static void printGrep(PsimagLite::String inputfile, ParametersForGrepType params)
	{
		SizeType           lenInput   = inputfile.size();
		PsimagLite::String dotcout    = ".cout";
		SizeType           lenDotcout = dotcout.size();
		SizeType           loc        = (lenInput < lenDotcout) ? 0 : lenInput - lenDotcout;
		bool               isCout     = (inputfile.substr(loc, lenDotcout) == dotcout);

		PsimagLite::String coutName
		    = (isCout) ? inputfile : ProgramGlobals::coutName(inputfile, "toolboxdmrg");
		std::ifstream fin(coutName.c_str());
		if (!fin || fin.bad()) {
			err("Could not open file " + coutName + "\n");
		}

		GrepForLabel::hook(fin, "", 1, params);
	}

	static void analize(const DmrgParametersType& solverParams,
	                    const GeometryType&       geometry,
	                    PsimagLite::String        extraOptions)
	{
		PairSizeStringType g = findLargestGeometry(geometry);
		SizeType           m = neededKeptStates(g, geometry, solverParams);
		std::cout << "Geometry= " << g << "\n";
		std::cout << "Needed m=" << m << "\n";
	}

private:

	static PairSizeStringType findLargestGeometry(const GeometryType& geometry)
	{
		SizeType terms = geometry.terms();
		assert(terms > 0);
		PsimagLite::String g            = geometry.term(0).label();
		SizeType           heaviestTerm = 0;
		for (SizeType i = 1; i < terms; ++i) {
			PsimagLite::String tmp = geometry.term(i).label();
			if (geometryGreater(tmp, g)) {
				g            = tmp;
				heaviestTerm = i;
			}
		}

		return PairSizeStringType(heaviestTerm, g);
	}

	static bool geometryGreater(PsimagLite::String g1, PsimagLite::String g2)
	{
		if (g1 == "longchain")
			return false;
		std::cerr << g1 << " " << g2 << "\n";
		return true;
	}

	static SizeType neededKeptStates(PairSizeStringType        g,
	                                 const GeometryType&       geometry,
	                                 const DmrgParametersType& solverParams)
	{
		SizeType m           = 0;
		SizeType modelFactor = getModelFactor(solverParams.model);
		SizeType n           = geometry.numberOfSites();
		if (g.second == "longchain") { // 1D
			return modelFactor * n; // modelFactor * Lx
		} else if (g.second == "ladder" || g.second == "ladderx") {
			SizeType Lx         = geometry.term(g.first).length(0);
			SizeType Ly         = geometry.term(g.first).length(1);
			SizeType TwoToTheLy = (1 << Ly);
			m = modelFactor * Lx * TwoToTheLy; // modelFactor * Lx * 2^Ly
		} else {
			err("neededKeptStates: unknown geometry" + g.second + "\n");
		}

		return m;
	}

	static SizeType getModelFactor(PsimagLite::String model)
	{
		if (model == "HubbardOneBand")
			return 7;
		if (model == "Heisenberg")
			return 4; // correct for s
		err("getModelFactor: unknown model" + model + "\n");
		return 0;
	}

}; // class ToolBox

} // namespace Dmrg
/*@}*/
#endif
