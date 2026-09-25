// SPDX-FileCopyrightText: Copyright (c) 2009-2016-2018-2022, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [DMRG++, Version 5.]
// [DMRG++ authors; see AUTHORS.]

/** \ingroup DMRG */
/*@{*/

/*! \file ProgramGlobals.h
 *
 *
 *
 */
#ifndef PROGRAM_LIMITS_H
#define PROGRAM_LIMITS_H
#include "../Version.h"
#include "Utils.h"
#include <PsimagLite/PsimagLite.h>
#include <PsimagLite/Vector.h>

#include <algorithm>
#include <numeric>

namespace Dmrg {

struct ProgramGlobals {

	static SizeType maxElectronsOneSpin;

	static bool oldChangeOfBasis;

	static const std::string license;

	static const SizeType MAX_LPS = 1000;

	enum class DirectionEnum
	{
		INFINITE,
		EXPAND_ENVIRON,
		EXPAND_SYSTEM
	};

	enum class ConnectionEnum
	{
		SYSTEM_SYSTEM,
		SYSTEM_ENVIRON,
		ENVIRON_SYSTEM,
		ENVIRON_ENVIRON
	};

	enum class SysOrEnvEnum
	{
		SYSTEM,
		ENVIRON
	};

	enum class FermionOrBosonEnum
	{
		FERMION,
		BOSON
	};

	enum class VerboseEnum
	{
		NO,
		YES
	};

	static FermionOrBosonEnum multipy(const FermionOrBosonEnum& a, const FermionOrBosonEnum& b)
	{
		if (a == FermionOrBosonEnum::BOSON)
			return b;

		return (b == FermionOrBosonEnum::BOSON) ? FermionOrBosonEnum::FERMION
		                                        : FermionOrBosonEnum::BOSON;
	}

	static void init(SizeType maxElectronsOneSpin_)
	{
		if (maxElectronsOneSpin == maxElectronsOneSpin_)
			return;
		if (maxElectronsOneSpin != 0) {
			std::cerr << PsimagLite::AnsiColor::blue;
			std::string msg("ProgramGlobals::init(...) replayed\n");
			std::cout << msg;
			std::cerr << msg;
			std::cerr << PsimagLite::AnsiColor::reset;
		}

		maxElectronsOneSpin = maxElectronsOneSpin_;
	}

	static int findBorderSiteFrom(SizeType site, DirectionEnum direction, SizeType n)
	{
		if (site == 1 && direction == DirectionEnum::EXPAND_ENVIRON)
			return 0;

		if (site == n - 2 && direction == DirectionEnum::EXPAND_SYSTEM)
			return n - 1;

		return -1;
	}

	static std::string rootName(std::string filename)
	{
		std::string rootname = filename;
		size_t      index    = rootname.find(".", 0);
		if (index != std::string::npos) {
			rootname.erase(index, filename.length());
		}

		return rootname;
	}

	static std::string coutName(std::string filename, std::string app_name)
	{
		std::string rootname = PsimagLite::basename(filename);
		size_t      index    = rootname.find(".", 0);
		if (index != std::string::npos) {
			rootname.erase(index, filename.length());
		}

		std::string maybe_app_name = (app_name == "DMRG++::observe") ? "observe" : "";

		return "runFor" + maybe_app_name + rootname + ".cout";
	}

	static SizeType logBase2(SizeType x)
	{
		SizeType counter = 0;
		while (x > 0) {
			x >>= 1;
			counter++;
		}

		return (counter == 0) ? counter : counter - 1;
	}

	static SizeType volumeOf(const PsimagLite::Vector<SizeType>::Type& v)
	{
		assert(v.size() > 0);
		SizeType ret = v[0];
		for (SizeType i = 1; i < v.size(); i++)
			ret *= v[i];
		return ret;
	}

	friend std::istream& operator>>(std::istream& is, DirectionEnum& direction)
	{
		int x = -1;
		is >> x;
		if (x == 0) {
			direction = DirectionEnum::INFINITE;
		} else if (x == 1) {
			direction = DirectionEnum::EXPAND_ENVIRON;
		} else if (x == 2) {
			direction = DirectionEnum::EXPAND_SYSTEM;
		} else {
			err("istream& operator>> DirectionEnum\n");
		}

		return is;
	}

	static std::string toString(const DirectionEnum d)
	{
		switch (d) {
		case DirectionEnum::INFINITE:
			return "INFINITE";
			break;
		case DirectionEnum::EXPAND_ENVIRON:
			return "EXPAND_ENVIRON";
			break;
		case DirectionEnum::EXPAND_SYSTEM:
			return "EXPAND_SYSTEM";
			break;
		}

		return "UNKNOWN_DIRECTION_ENUM";
	}

	static std::string killSpaces(std::string str)
	{
		std::string    buffer;
		const SizeType n = str.length();
		for (SizeType i = 0; i < n; ++i)
			if (str[i] != ' ')
				buffer += str[i];
		return buffer;
	}

	static std::string toLower(std::string data)
	{
		std::transform(data.begin(),
		               data.end(),
		               data.begin(),
		               [](unsigned char c) { return std::tolower(c); });
		return data;
	}

	static std::string SYSTEM_STACK_STRING;
	static std::string ENVIRON_STACK_STRING;
}; // ProgramGlobals

} // namespace Dmrg
/*@}*/
#endif
