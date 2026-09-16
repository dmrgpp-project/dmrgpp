// SPDX-FileCopyrightText: Copyright (c) 2009-2017, UT-Battelle, LLC
// SPDX-License-Identifier: GPL-3.0-or-later
//
// [PsimagLite, Version 1.]
// [by G.A., Oak Ridge National Laboratory]

/** \ingroup PsimagLite */
/*@{*/

/*! \file HostInfo.h
 *
 * Information about the host computer
 */

#ifndef APPLICATION_INFO_H
#define APPLICATION_INFO_H

#include "BitManip.h"
#include "Io/IoSerializerStub.h"
#include "MersenneTwister.h"
#include "Vector.h"
#include <cassert>
#include <cstdlib>
#include <sys/time.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

namespace PsimagLite {

class ApplicationInfo {

public:

	using RunIdType = String;

	ApplicationInfo(const PsimagLite::String& name)
	    : name_(name)
	    , pid_(getpid())
	    , runId_(runIdInternal())
	    , isFinalized_(false)
	{ }

	void finalize() { isFinalized_ = true; }

	time_t unixTime(bool arg = false) const
	{
		struct timeval tv;
		gettimeofday(&tv, 0);
		return (arg) ? tv.tv_usec : tv.tv_sec;
	}

	String getTimeDate() const
	{
		time_t tt = unixTime();
		return asctime(localtime(&tt));
	}

	String hostname() const
	{
		int    len  = 1024;
		char*  name = new char[len];
		int    ret  = gethostname(name, len);
		String retString;
		if (ret != 0) {
			retString = "UNKNOWN";
		} else {
			retString = name;
		}

		delete[] name;

		return retString;
	}

	const RunIdType runId() const { return runId_; }

	unsigned int pid() const { return pid_; }

	void write(String label, IoSerializer& serializer) const
	{
		String root = label;
		if (!isFinalized_) {
			serializer.createGroup(root);

			serializer.write(root + "/Name", name_);
			serializer.write(root + "/RunId", runId_);
			serializer.write(root + "/UnixTimeStart", unixTime(false));
		} else {
			serializer.write(root + "/UnixTimeEnd", unixTime(false));
		}
	}

	static void setEnv(String name, String value)
	{
		int ret = setenv(name.c_str(), value.c_str(), true);
		if (ret != 0)
			throw RuntimeError("Could not setenv " + name + "=" + value + "\n");
		std::cout << "Set " << name << "=" << value << "\n";
	}

	friend std::ostream& operator<<(std::ostream& os, const ApplicationInfo& ai)
	{
		if (ai.isFinalized_)
			printFinalLegacy(os, ai);
		else
			printInit(os, ai);

		return os;
	}

private:

	static void printInit(std::ostream& os, const ApplicationInfo& ai)
	{
		os << ai.getTimeDate();
		os << "Hostname: " << ai.hostname() << "\n";
		os << "RunID=" << ai.runId_ << "\n";
		os << "UnixTimeStart=" << ai.unixTime(false) << "\n";
		os << "SizeType=" << sizeof(SizeType) << "\n";
	}

	static void printFinalLegacy(std::ostream& os, const ApplicationInfo& ai)
	{
		OstringStream msg(std::cout.precision());
		msg() << ai.name_ << "\nsizeof(SizeType)=" << sizeof(SizeType) << "\n";
#ifdef USE_FLOAT
		msg() << ai.name_ << " using float\n";
#else
		msg() << ai.name_ << " using double\n";
#endif
		msg() << "UnixTimeEnd=" << ai.unixTime(false) << "\n";
		msg() << ai.getTimeDate();
		os << msg().str();
	}

	RunIdType runIdInternal() const
	{
		unsigned int                      p  = getpid();
		time_t                            tt = unixTime(true);
		MersenneTwister                   mt(tt + p);
		unsigned int                      x = tt ^ mt.random();
		OstringStream                     msgg(std::cout.precision());
		OstringStream::OstringStreamType& msg = msgg();
		msg << x;
		x = p ^ mt.random();
		msg << x;
		unsigned long int y = atol(msg.str().c_str());
		y ^= mt.random();
		x = BitManip::countKernighan(y);
		OstringStream                     msgg2(std::cout.precision());
		OstringStream::OstringStreamType& msg2 = msgg2();
		msg2 << y;
		if (x < 10)
			msg2 << "0";
		msg2 << x;
		return msg2.str();
	}

	long unsigned int convertToLuint(PsimagLite::String str) const
	{
		long unsigned int sum  = 0;
		long unsigned int prod = 1;
		int               l    = str.length();
		assert(l < 20);

		for (int i = 0; i < l; ++i) {
			unsigned int c = str[l - i - 1] - 48;
			sum += prod * c;
			prod *= 10;
		}

		return sum;
	}

	PsimagLite::String name_;
	unsigned int       pid_;
	const RunIdType    runId_;
	bool               isFinalized_;
}; // class ApplicationInfo

std::ostream& operator<<(std::ostream& os, const ApplicationInfo& ai);

} // namespace PsimagLite

/*@}*/
#endif // APPLICATION_INFO_H
